!===================
  module data_posteriori
!===================

  use global_parameters, ONLY : realkind

  implicit none

  public 

! include 'mesh_params.h'
! include 'mesh_params_kernel.h'

  integer :: num_model_param,num_model_reconst,num_windows,need_filt
  character(len=11), allocatable :: filter_char(:)
  character(len=15),allocatable :: src_type1fwd_post(:),src_type2fwd_post(:)
  character(len=15),allocatable :: src_type1bwd_post(:),src_type2bwd_post(:)
  character(len=10),allocatable :: fwdstf_type_post(:),bwdstf_type_post(:)
  real(kind=realkind) :: Mij_post(6)
  logical :: xc_tt,xc_am,pk_tt,on_tt,insta,wf_dd,phase,envel,model_norm,do_filt
  logical :: need_omkern,hilbert_needed,gabor_needed,instph,insten
  logical :: do_lam_post,do_mu_post,do_rho_post,do_vp_post,do_vs_post,do_imped_post
  logical :: do_vbulk_post,do_vsbulk_post,do_rhobulk_post
  integer, allocatable :: dim_fwd_post(:),dim_bwd_post(:),dump_ind(:,:),pick_ind(:,:)
  integer, allocatable :: num_picks(:)
  logical, allocatable :: travelpick(:)
  real(kind=realkind), allocatable :: post_per(:,:),dump_time(:,:),time_pick(:,:),time_pick_read(:)
  real(kind=realkind), allocatable :: Mijsum(:),misfit_prefactor(:,:,:),norm_fac(:,:,:)
  real(kind=realkind), allocatable :: ttdelay(:,:),ampdiff(:,:)
  real(kind=realkind), allocatable :: wft1(:),wft2(:),a0(:)
  character(len=5),allocatable :: misfit(:),pick_misfit(:)
  character(len=4),allocatable :: filter_type_posteriori(:)
  character(len=3),allocatable :: model_run(:),model_reconst(:),wfkernel_type(:)
  character(len=3),allocatable :: app_period1(:),app_period2(:)
  integer, allocatable :: reccomp_post(:)
  real(kind=realkind), allocatable :: u0post(:),v0post(:)
  
  !=======================
  end module data_posteriori
!=======================
