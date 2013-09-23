!====================
module coord_trafo
!====================

  !use data_mesh
  use global_parameters
  !use data_rota

  implicit none

  !public :: init_rotations,get_r_theta,get_r_theta_dble
  public :: xyz2rthetaphi,xyz2rthetaphi_dble,rthetaphi2xyz !,rotate_velocity_components
  public :: sphi2xy_arr,sphi2xy,xyz2sphiz
!  public :: azimuthal_prefactor_fwd_isim
 ! public :: rotate_strain_components,rotate_vec_cyl2cart,rotate_tens_cyl2cart
  private

  contains 
!/////////////////////////////////////////////////
!
!!-------------------------------------------------------------------------
!subroutine init_rotations(thr,phr)
!  real(kind=realkind)          :: thr,phr
!
!  if (mynum==0)write(6,*)mynum,' do some rotations........'; call flush(6)
!  call azimuthal_prefactor
!
!  if (mynum==0)write(6,*)mynum,' compute receiver rotation matrix........'; call flush(6)
!  call compute_rec_rotmatrix(thr,phr)
!
!end subroutine init_rotations
!!-------------------------------------------------------------------------
!
!!-------------------------------------------------------------------------
!subroutine compute_rec_rotmatrix(thr,phir)
!
!  include 'mesh_params.h'
!  include 'mesh_params_kernel.h'
!
!  real(kind=realkind), intent(in) :: thr,phir
!  integer :: i,j
!
!! This is the rotation matrix to transform tensors from
!! rotated xyz to the unrotated xyz
!! column 1
!    vec_rot(1,1) = cos(thr)*cos(phir)
!    vec_rot(2,1) = cos(thr)*sin(phir)
!    vec_rot(3,1) = -sin(thr)
!
!! column2
!    vec_rot(1,2) = -sin(phir)
!    vec_rot(2,2) = cos(phir)
!    vec_rot(3,2) = 0.d0
!
!! column3
!    vec_rot(1,3) = sin(thr)*cos(phir)
!    vec_rot(2,3) = sin(thr)*sin(phir)
!    vec_rot(3,3) = cos(thr)
!
!! transforms from unrotated to rotated xyz
!    trans_vec_rot=transpose(vec_rot)
!
!    write(6,*)'vector rotation matrix:'
!    do i=1,3
!       write(6,*)(vec_rot(i,j),j=1,3)
!    enddo
!
!end subroutine compute_rec_rotmatrix
!!-------------------------------------------------------------------------
!
!!-------------------------------------------------------------------------
!subroutine azimuthal_prefactor
!
!  use data_arrays, only : do_rho,do_lam,do_mu
!  use input_output, only : save_kernel
!
!  character(len=200) :: fname
!
!  allocate(azim1_fwd(npts),azim2_fwd(npts),azim1_bwd(npts),azim2_bwd(npts))
!
!  if (src_type1fwd/='monopole' .or. src_type1bwd/='monopole') then
!     do_azim=.true.
!  else
!     do_azim=.false.
!     if (mynum==0) write(6,*)"don't need any azimuthal prefactor: fwd/bwd monopoles"
!  endif
!
!! defaults
!     azim1_fwd = 1.d0; if (allocated(azim2_fwd)) azim2_fwd = 0.d0
!     azim1_bwd = 1.d0; if (allocated(azim2_bwd)) azim2_bwd = 0.d0
!
!     if (nsim==1) call azimuthal_prefactor_fwd_isim(1,azim1_fwd,azim2_fwd)
!
!     ! backward fields
!     if (src_type2bwd=='mxz' .or. src_type2bwd=='xforce') then
!        if (mynum==0) write(6,*) &
!             '  Calculating prefactors for correct azimuth for',src_type2bwd
!        azim1_bwd =  cos(phird)
!        azim2_bwd =  -sin(phird)
!
!     elseif (src_type2bwd=='myz' .or. src_type2bwd=='yforce') then 
!        if (mynum==0) write(6,*) &
!             '  Calculating prefactors for correct azimuth for',src_type2bwd
!        azim1_bwd = sin(phird)
!        azim2_bwd=  cos(phird)
!        
!     elseif (src_type2bwd=='mxx_m_myy') then 
!        if (mynum==0) write(6,*) &
!             '  Calculating prefactors for correct azimuth for',src_type2bwd
!        azim1_bwd =  cos(2.d0*phird)
!        azim2_bwd  =  -sin(2.d0*phird)
!        
!     elseif (src_type2bwd=='mxy') then 
!        if (mynum==0) write(6,*) &
!             '  Calculating prefactors for correct azimuth for',src_type2bwd
!        azim1_bwd =  sin(2.d0*phird)
!        azim2_bwd =  cos(2.d0*phird)
!     endif
!
!     allocate(sinph(1:npts),cosph(1:npts))
!     sinph = sin(phird)
!     cosph = cos(phird)
!
!! plot phird and phigd
!     phird=phird*180./pi; phigd=phigd*180./pi
!     fname='phird_'; call save_kernel(npts,phird,fname,appmynum)
!     fname='phigd_'; call save_kernel(npts,phigd,fname,appmynum)
!     fname='az1fwd'; call save_kernel(npts,azim1_fwd,fname,appmynum)
!     fname='az2fwd'; call save_kernel(npts,azim2_fwd,fname,appmynum)
!     fname='az1bwd'; call save_kernel(npts,azim1_bwd,fname,appmynum)
!     fname='az2bwd'; call save_kernel(npts,azim2_bwd,fname,appmynum)
!
!     deallocate(phird,phigd)
!
!end subroutine azimuthal_prefactor
!!-------------------------------------------------------------------------
!
!!-----------------------------------------------------------------------------------------
!subroutine azimuthal_prefactor_fwd_isim(isim,az1,az2)
!  implicit none
!  integer, intent(in) :: isim
!  integer :: i
!  real(kind=realkind),intent(out) :: az1(1:npts),az2(1:npts)
!
!  write(6,*)'azimuthal prefactor fwd isim:',isim
!
!! forward fields in frame gd
!  select case(isim)
!  case(1) ! Mzz
!     az1(1:npts) = Mij(1)
!     az2 = 0.d0
!     
!  case(2) ! Mxx+Myy
!     az1=0.5*(Mij(2)+Mij(3))
!     az2 = 0.d0
!     
!  case(3) ! dipole
!     az1 = Mij(4)*cos(phigd) + Mij(5)*sin(phigd)
!     az2 = Mij(5)*cos(phigd) - Mij(4)*sin(phigd) 
!     
!  case(4) ! quadrupole
!     az1 = 0.5*(Mij(2)-Mij(3))*cos(2.d0*phigd) + Mij(6)*sin(2.d0*phigd)
!     az2 = Mij(6)*cos(2.d0*phigd) - 0.5*( Mij(2)-Mij(3) )*sin(2.d0*phigd) 
!
!  case default
!     write(6,*)mynum,'unknown number of simulations',isim
!  end select
!  
!  write(6,*)'done with azimuthal prefactor fwd isim.'
!
!end subroutine azimuthal_prefactor_fwd_isim
!!-----------------------------------------------------------------------------------------
!
!!-------------------------------------------------------------------------
!subroutine rotate_vec_cyl2cart(uin,uout,n,azim1,azim2)
!
!use data_rota
!
!integer, intent(in) :: n
!real(kind=realkind), intent(in) :: uin(npts,n)
!real(kind=realkind), intent(out) :: uout(npts,3)
!integer :: i
!real(kind=realkind) :: azim1(npts),azim2(npts)
!
!write(6,*)'AZIM:',maxval(azim1),maxval(azim2),maxval(phird)
!
!if (n==3) then 
!   uout(:,1) = azim1*cosph*uin(:,1) - azim2*sinph*uin(:,2)
!   uout(:,2) = azim1*sinph*uin(:,1) + azim2*cosph*uin(:,2)
!else
!   uout(:,1)= azim1*cosph*uin(:,1)
!   uout(:,2) = 0.
!endif
!uout(:,3)=azim1*uin(:,n)
!
!end subroutine rotate_vec_cyl2cart
!!-------------------------------------------------------------------------
!
!!-------------------------------------------------------------------------
!subroutine rotate_velocity_components(uout)
!
!use data_rota
!
!real(kind=realkind), intent(inout) :: uout(1:npts,1:ndim)
!integer :: ipt
!real(kind=realkind) :: uouttmp(3)
!
!do ipt=1,npts
!   uouttmp=uout(ipt,1:ndim)
!   uout(ipt,1:ndim)=matmul(vec_rot,uouttmp)
!enddo
!
!end subroutine rotate_velocity_components
!!-------------------------------------------------------------------------
!
!!-----------------------------------------------------------------------------------------
!subroutine rotate_tens_cyl2cart(nnpts,uin,uout,azim1,azim2,n)
!
!use data_rota
!
!integer, intent(in) :: n,nnpts
!real(kind=realkind), intent(in) :: uin(nnpts,n),azim1(nnpts),azim2(nnpts)
!real(kind=realkind), intent(out) :: uout(nnpts,6)
!integer :: i
!real(kind=realkind) :: matr(3,3),field(3,3),fieldout(3,3),transmatr(3,3)
!
!do i=1,nnpts
!   field=0.
!   field(1,1)=uin(i,1);                      field(1,3)=uin(i,4); ! Esz
!                        field(2,2)=uin(i,2)  
!   field(3,1)=uin(i,4);                      field(3,3)=uin(i,3);
!
!   if (n==6) then 
!      field(1,2)=uin(i,5); field(2,1)=uin(i,5) !Esphi
!      field(2,3)=uin(i,6); field(3,2)=uin(i,6) !Ezphi
!   endif
!
!   ! hadamard product: azimuthal dependence
!   matr(1,1)=azim1(i) ; matr(1,2)=azim2(i) ;matr(1,3)=azim1(i)
!   matr(2,1)=azim2(i) ; matr(2,2)=azim1(i) ;matr(2,3)=azim2(i)
!   matr(3,1)=azim1(i) ; matr(3,2)=azim2(i) ;matr(3,3)=azim1(i) 
!   field=matr*field
!
!   ! rotation to xyz in receiver frame
!   matr=0.
!   matr(1,1)=cosph(i); matr(1,2)=sinph(i); 
!   matr(2,1)=-sinph(i); matr(2,2)=cosph(i); matr(3,3)=1.
!   transmatr=transpose(matr); fieldout=matmul(matr,field); fieldout=matmul(fieldout,transmatr)
!
!! symmetric strain?
!   if (abs(fieldout(1,3)-fieldout(3,1))/abs(fieldout(1,3)) >1.e-4) then 
!      write(6,*)'nonsymmetric strain components (1,3) at',i
!      write(6,*)'(1,3),(3,1):',fieldout(1,3),fieldout(3,1)
!!      stop
!   endif
!   if (n==6) then
!      if (abs(fieldout(2,3)-fieldout(3,2))/abs(fieldout(2,3)) >1.e-4) then 
!         write(6,*)'nonsymmetric strain components (2,3) at',i
!         write(6,*)'(2,3),(3,2):',fieldout(2,3),fieldout(3,2)
!!         stop
!      elseif (abs(fieldout(1,2)-fieldout(2,1))/abs(fieldout(1,2))>1.e-4) then 
!         write(6,*)'nonsymmetric strain components (1,2) at',i
!         write(6,*)'(1,2),(2,1):',fieldout(1,2),fieldout(2,1)
!!         stop
!      endif
!   endif
!   
!   ! back to vector
!   uout(i,1)=fieldout(1,1); uout(i,2)=fieldout(2,2)
!   uout(i,3)=fieldout(3,3); uout(i,4)=fieldout(1,3)
!   if (n==6) then 
!      uout(i,5)=fieldout(1,2); uout(i,6)=fieldout(2,3)
!   endif
!enddo
!
!end subroutine rotate_tens_cyl2cart
!!-------------------------------------------------------------------------
!
!!-------------------------------------------------------------------------
!subroutine rotate_strain_components(nnpts,uin)
!
!use data_rota
!integer, intent(in) :: nnpts
!real(kind=realkind), intent(inout) :: uin(nnpts,6)
!integer :: i
!real(kind=realkind) :: matr(3,3),field(3,3),fieldout(3,3)
!real :: eps,maxfield
!
!eps=1.e-3
!
!if (ndim/=6) then
!   write(6,*)'issue... rotating strain components of field of size:',ndim
!   stop
!endif
!
!
!do i=1,nnpts
!   field=0.; fieldout=0.
!   field(1,1)=uin(i,1);  field(2,1)=uin(i,5); field(3,1)=uin(i,4)
!   field(1,2)=uin(i,5);  field(2,2)=uin(i,2); field(3,2)=uin(i,6)
!   field(1,3)=uin(i,4);  field(2,3)=uin(i,6); field(3,3)=uin(i,3)
!
!   ! rotation to global xyz
!   fieldout=matmul(vec_rot,field); field=matmul(fieldout,trans_vec_rot)
!
!! symmetric strain?
!   maxfield = max(1.e-15,0.01*maxval(abs(field)))
!   if ( abs(field(1,3)) > maxfield .and. abs(field(1,3)-field(3,1))>abs(field(1,3))*eps) then 
!      write(6,*)'nonsymmetric strain components (1,3) at',i
!      write(6,*)'(1,3),(3,1):',field(1,3),field(3,1),maxfield
!      stop
!   elseif (abs(field(1,3)) > maxfield .and. abs(field(2,3)-field(3,2))>abs(field(2,3))*eps) then 
!      write(6,*)'nonsymmetric strain components (2,3) at',i
!      write(6,*)'(2,3),(3,2):',field(2,3),field(3,2),maxfield
!      stop
!   elseif (abs(field(1,3)) > maxfield .and. abs(field(1,2)-field(2,1))>abs(field(1,2))*eps) then 
!      write(6,*)'nonsymmetric strain components (1,2) at',i
!      write(6,*)'(1,2),(2,1):',field(1,2),field(2,1),maxfield
!      stop
!   endif
!
!   ! back to vector
!   uin(i,1)=field(1,1) ;  uin(i,2)=field(2,2);   uin(i,3)=field(3,3)
!   uin(i,4)=field(3,1);   uin(i,5)=field(2,1);   uin(i,6)=field(3,2)
!enddo
!
!
!end subroutine rotate_strain_components
!!-------------------------------------------------------------------------
!
!!-------------------------------------------------------------------------
!subroutine matmul_test
!implicit none
!
!integer :: i,j,k,l,in,im
!real(kind=realkind) :: rot(3,3),e(3,3),erot(3,3),tmp(3,3),tens_rot(3,3),rott(3,3)
!real(kind=realkind) :: rot2(3,3),rot2t(3,3),tens_rot2(3,3),A(3,3),joint_mat(6,6)
!
!real(kind=realkind) :: tens6_new(6),tens6_old(6)
!
!!rot(1,1)= 2; rot(1,2)=3; rot(1,3)=0;
!!rot(2,1)= -3; rot(2,2)=2; rot(2,3)=0;
!!rot(3,1)= 0; rot(3,2)=0; rot(3,3)=1;
!
!!rot(1,1)= 0; rot(1,2)=1; rot(1,3)=0;
!!rot(2,1)= -1; rot(2,2)=0; rot(2,3)=0;
!!rot(3,1)= 0; rot(3,2)=0; rot(3,3)=1;
!
!rot(1,1)= 1; rot(1,2)=0; rot(1,3)=0;
!rot(2,1)= 0; rot(2,2)=1; rot(2,3)=0;
!rot(3,1)= 0; rot(3,2)=0; rot(3,3)=1;
!
!e(:,:)=0; erot(:,:)=0; tmp(:,:)=0
!e(1,1)= 10; e(1,2)= 100; e(1,3)= 1000;
!e(2,2)= 10000; e(3,3)= 100000; e(2,3)= 1000000;
!e(2,1)=e(1,2); e(3,1)=e(1,3); e(3,2)=e(2,3); 
!
!rot2(1,1)= 0.1111; rot2(1,2)=0.2222; rot2(1,3)=0.3333;
!rot2(2,1)= 0.4444; rot2(2,2)=0.5555; rot2(2,3)=0;
!rot2(3,1)= 0.6666; rot2(3,2)=0.7777; rot2(3,3)=0.8888;
!
!rott=transpose(rot);
!rot2t=transpose(rot2);
!tmp=matmul(rott,e);
!erot=matmul(tmp,rot);
!
!write(6,*)'rotation matrix:'
!write(6,*)rot(1,1), rot(1,2), rot(1,3)
!write(6,*)rot(2,1), rot(2,2), rot(2,3)
!write(6,*)rot(3,1), rot(3,2), rot(3,3)
!
!write(6,*)'input symmetric matrix e:'
!write(6,*)e(1,1), e(1,2), e(1,3)
!write(6,*)e(2,1), e(2,2), e(2,3)
!write(6,*)e(3,1), e(3,2), e(3,3)
!
!write(6,*)
!write(6,*)'rot^T  e:'
!write(6,*)tmp(1,1), tmp(1,2), tmp(1,3)
!write(6,*)tmp(2,1), tmp(2,2), tmp(2,3)
!write(6,*)tmp(3,1), tmp(3,2), tmp(3,3)
!
!write(6,*)
!write(6,*)'output matrix (rot^T e rot):'
!write(6,*)erot(1,1), erot(1,2), erot(1,3)
!write(6,*)erot(2,1), erot(2,2), erot(2,3)
!write(6,*)erot(3,1), erot(3,2), erot(3,3)
!
!! second rotation
!tmp=matmul(rot2t,erot);
!erot=matmul(tmp,rot2);
!
!write(6,*)
!write(6,*)'R^T (rot^T e rot) R:'
!write(6,*)erot(1,1), erot(1,2), erot(1,3)
!write(6,*)erot(2,1), erot(2,2), erot(2,3)
!write(6,*)erot(3,1), erot(3,2), erot(3,3)
!
!! check commutativity:
!
!tmp=matmul(rot2t,e);
!erot=matmul(tmp,rot2);
!tmp=matmul(rott,erot);
!erot=matmul(tmp,rot);
!
!write(6,*)
!write(6,*)'rot^T (R^T  e  R) rot):'
!write(6,*)erot(1,1), erot(1,2), erot(1,3)
!write(6,*)erot(2,1), erot(2,2), erot(2,3)
!write(6,*)erot(3,1), erot(3,2), erot(3,3)
!
!tmp=matmul(rot2t,rot2);
!erot=matmul(tmp,e);
!tmp=matmul(rott,rot);
!erot=matmul(tmp,erot);
!
!write(6,*)
!write(6,*)'rot^T rot (R^T R e  R) ):'
!write(6,*)erot(1,1), erot(1,2), erot(1,3)
!write(6,*)erot(2,1), erot(2,2), erot(2,3)
!write(6,*)erot(3,1), erot(3,2), erot(3,3)
!
!tens_rot=0
!! alternatively:
!   do i=1,3
!      do j=1,3
!!      diagonal contributions
!         do k=1,3
!            tens_rot(i,j) = tens_rot(i,j) + rott(i,k)*rott(j,k)*e(k,k)
!          enddo
!!      off-diagonal (symmtric!) contributions
!          tens_rot(i,j) = tens_rot(i,j) + &
!                          (rott(i,1)*rott(j,2) + rott(i,2)*rott(j,1))*e(1,2)+&
!                          (rott(i,1)*rott(j,3) + rott(i,3)*rott(j,1))*e(1,3)+&
!                          (rott(i,2)*rott(j,3) + rott(i,3)*rott(j,2))*e(2,3)
!      enddo
!   enddo
!
!write(6,*)
!write(6,*)'output matrix (rot^T e rot) unrolled:'
!write(6,*)tens_rot(1,1), tens_rot(1,2), tens_rot(1,3)
!write(6,*)tens_rot(2,1), tens_rot(2,2), tens_rot(2,3)
!write(6,*)tens_rot(3,1), tens_rot(3,2), tens_rot(3,3)
!
!tens_rot2=0.0
!! second rotation
!   do i=1,3
!      do j=1,3
!!      diagonal contributions
!         do k=1,3
!            tens_rot2(i,j) = tens_rot2(i,j) + rot2t(i,k)*rot2t(j,k)*tens_rot(k,k)
!          enddo
!!      off-diagonal (symmtric!) contributions
!          tens_rot2(i,j) = tens_rot2(i,j) + &
!                (rot2t(i,1)*rot2t(j,2) + rot2t(i,2)*rot2t(j,1))*tens_rot(1,2)+&
!                (rot2t(i,1)*rot2t(j,3) + rot2t(i,3)*rot2t(j,1))*tens_rot(1,3)+&
!                (rot2t(i,2)*rot2t(j,3) + rot2t(i,3)*rot2t(j,2))*tens_rot(2,3)
!      enddo
!   enddo
!
!write(6,*)
!write(6,*)'output matrix R^T(rot^T e rot)R unrolled:'
!write(6,*)tens_rot2(1,1), tens_rot2(1,2), tens_rot2(1,3)
!write(6,*)tens_rot2(2,1), tens_rot2(2,2), tens_rot2(2,3)
!write(6,*)tens_rot2(3,1), tens_rot2(3,2), tens_rot2(3,3)
!
!! check condensed version
!A = 0
!do i=1,3
!   do j=1,3
!      do k=1,3
!         A(i,j) = A(i,j) + rott(i,k) * rot2t(k,j)
!      enddo
!   enddo
!enddo
!
!tens_rot2=0
!   do i=1,3
!      do j=1,3
!!      diagonal contributions
!         do k=1,3
!            tens_rot2(i,j) = tens_rot2(i,j) + A(i,k)*A(j,k)*e(k,k)
!          enddo
!!      off-diagonal (symmtric!) contributions
!          tens_rot2(i,j) = tens_rot2(i,j) + &
!                (A(i,1)*A(j,2) + A(i,2)*A(j,1))*e(1,2)+&
!                (A(i,1)*A(j,3) + A(i,3)*A(j,1))*e(1,3)+&
!                (A(i,2)*A(j,3) + A(i,3)*A(j,2))*e(2,3)
!      enddo
!   enddo
!
!write(6,*)
!write(6,*)'output matrix A^T e A, A=rot^T R^T unrolled:'
!write(6,*)tens_rot2(1,1), tens_rot2(1,2), tens_rot2(1,3)
!write(6,*)tens_rot2(2,1), tens_rot2(2,2), tens_rot2(2,3)
!write(6,*)tens_rot2(3,1), tens_rot2(3,2), tens_rot2(3,3)
!
!
!! do it jointly via 
!joint_mat = 0.d0
!    do in=1,6
!       if (in<=3) then 
!          i=in; j=in
!       elseif (in==4) then 
!          i=1; j=2    
!       elseif (in==5) then 
!          i=1; j=3
!       elseif (in==6) then 
!          i=2; j=3
!       endif         
!       do im=1,3
!          joint_mat(in,im) = A(i,im)*A(j,im) ! diagonal elements
!       enddo
!       joint_mat(in,4) = A(i,1)*A(j,2) + A(i,2)*A(j,1) !12
!       joint_mat(in,5) = A(i,1)*A(j,3) + A(i,3)*A(j,1) !13
!       joint_mat(in,6) = A(i,2)*A(j,3) + A(i,3)*A(j,2) !23
!    enddo
!
!do i=1,3
!   tens6_old(i) = e(i,i)
!enddo
!tens6_old(4) = e(1,2); tens6_old(5) = e(1,3); tens6_old(6) = e(2,3); 
!
!do i=1,6
!   do j=1,6
!      tens6_new(i) = tens6_new(i) + joint_mat(i,j)* tens6_old(j)
!   enddo
!enddo
!
!write(6,*)
!write(6,*)'output matrix E(i) = joint_mat(ij) eold(j) :'
!write(6,*)tens6_new(1), tens6_new(4), tens6_new(5)
!write(6,*)tens6_new(4),tens6_new(2),tens6_new(6)
!write(6,*)tens6_new(5),tens6_new(6),tens6_new(3)
!
!end subroutine matmul_test
!!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
 subroutine get_r_theta(s,z,r,th)

real(kind=realkind), intent(in) :: s,z
real(kind=realkind), intent(out) :: r,th

! th=atan(s/(z+epsi))
th=atan2(s,z)

 if ( 0.0 > th ) th = real(pi) + th
 if (th == zero .and. z < 0) th = real(pi)

 r=sqrt(s**2 + z**2)

end subroutine get_r_theta
!==========================================================================
!

!-------------------------------------------------------------------------
 subroutine get_r_theta_dble(s,z,r,th)

double precision, intent(in) :: s,z
double precision, intent(out) :: r,th

 th=datan(s/(z+epsi))
 if ( 0.d0 > th ) th = dble(pi) + th
 if (th == zero .and. z < 0.d0) th = dble(pi)
 r=dsqrt(s**2 + z**2)

end subroutine get_r_theta_dble
!==========================================================================


!dk compute_slphilzl-------------------------------------------------------
  subroutine compute_slphilzl(sl,phil,zl,xgd,ygd,zgd,xorig,yorig,zorig)
  real(kind=realkind), intent(out)::  sl,phil,zl
  real(kind=realkind), intent(in) :: xgd,ygd,zgd
  real(kind=realkind), intent(in) :: xorig,yorig,zorig
  real(kind=realkind) :: arr
! S =(1/rr) * [ (yd*zr-zd*yr)^2 + (zd*xr-xd*zr)^2 + (xd*yr-yd*xr)^2 ]^(1/2)
! and
! Z = (1/rr) * (xd*xr+yd*yr+zd*zr)
  arr = (xorig**2+yorig**2+zorig**2)**(-.5)
!
  sl  =arr*( (ygd*zorig-zgd*yorig)**2 + &
             (zgd*xorig-xgd*zorig)**2 + &
             (xgd*yorig-ygd*xorig)**2     )**(.5)
  phil=0.
  zl = arr*(xgd*xorig+ygd*yorig+zgd*zorig)
!
end subroutine compute_slphilzl
!--------------------------------------------------------------------------




!dk rthetaphi2xyz----------------------------------------------------------
  subroutine rthetaphi2xyz(x,y,z,r,theta,phi)
  real(kind=realkind), intent(out) :: x,y,z
  real(kind=realkind), intent(in) :: r,theta,phi
  x=r*sin(theta)*cos(phi)
  y=r*sin(theta)*sin(phi)
  z=r*cos(theta)
  end subroutine rthetaphi2xyz
!--------------------------------------------------------------------------

!dk szphi2xyz----------------------------------------------------------
  subroutine sphi2xy(x,y,s,phi)
  real(kind=realkind), intent(out) :: x,y
  real(kind=realkind), intent(in) :: s,phi

  x=s*cos(phi)
  y=s*sin(phi)

  end subroutine sphi2xy
!--------------------------------------------------------------------------

!dk szphi2xyz_arr----------------------------------------------------------
  subroutine sphi2xy_arr(x,y,s,phi,n)
    integer, intent(in) :: n
    real(kind=realkind), intent(out) :: x(n),y(n)
    real(kind=realkind), intent(in) :: s(n),phi
    integer :: i

    do i=1,n
       x(i)=s(i)*cos(phi)
       y(i)=s(i)*sin(phi)
    enddo
    
  end subroutine sphi2xy_arr
!--------------------------------------------------------------------------



!dk xyz2rthetaphi----------------------------------------------------------
  subroutine xyz2rthetaphi(r,theta,phi,x,y,z)
  real(kind=realkind), intent(out) :: r,theta,phi
  real(kind=realkind), intent(in) :: x,y,z
  real(kind=realkind) :: pi
  pi = 2.*asin(1.)
  r = sqrt(x**2+y**2+z**2)
  theta = .5*pi-asin(z/r)

  phi = atan2(y,x)
  if (phi<0.0) phi=phi+2.d0*pi

  end subroutine xyz2rthetaphi
!--------------------------------------------------------------------------
!

!dk xyz2rthetaphi----------------------------------------------------------
  subroutine xyz2rthetaphi_dble(r,theta,phi,x,y,z)
  double precision, intent(out) :: r,theta,phi
  double precision, intent(in) :: x,y,z
  double precision :: pi
  pi = 2.d0*dasin(1.d0)
  r = dsqrt(x**2+y**2+z**2)
  theta = .5d0*pi-dasin(z/r)
  if (y>0) then
   if (x>0) then
    phi = datan(y/(x+1.e-20))
   else
    phi = pi+datan(y/(x+1.e-20))
   endif
  else
   if (x>0) then
    phi = 2.d0*pi+datan(y/(x+1.e-20))
   else
    phi = pi+datan(y/(x+1.e-20))
   end if
  end if
  if (dabs(x)<1.e-20) then
   if(y>0.) then
    phi = .5*pi
   else
    phi = 1.5*pi
   end if
  endif
  end subroutine xyz2rthetaphi_dble
!--------------------------------------------------------------------------

!dk xyz2sphiz--------------------------------------------------------------
  subroutine xyz2sphiz(s,phi,x,y,z)
  real(kind=realkind), intent(out) :: s, phi
  real(kind=realkind), intent(in) :: x,y,z
  real(kind=realkind) :: pi
  pi = 2.*asin(1.)
  s  = sqrt(x**2+y**2)

phi=atan2(y,x)
if (phi<0.d0) phi=2.d0*pi+phi

  end subroutine xyz2sphiz
!--------------------------------------------------------------------------
!


!====================
end module coord_trafo
!====================
