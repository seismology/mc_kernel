!===================
module data_mesh
!===================

  use global_parameters, ONLY : realkind

  implicit none

  public 

  include 'mesh_params.h'
  include 'mesh_params_kernel.h'

  integer                                        :: mynum, nproc
  character(len=4)                               :: appmynum, appnproc, char_wfkern, td_fd
  integer                                        :: ierror, reccomp
  integer                                        :: dim_fwd, dim_bwd, ndim
  integer                                        :: nsize, npts, n_src, nloc 
  integer                                        :: nsize_flu, nsize_sol 
  real(kind=realkind), dimension(ndumps)         :: time, timestep
  integer, dimension(nelem)                      :: ielr
  integer, dimension(nelem,2)                    :: proc_pt
  logical                                        :: cell_topology 
 
  real                                           :: globmax_lam, globmax_mu, globmax_rho
  real                                           :: globmax_vp, globmax_vs, globmax_imp

  character(len=100),allocatable                 :: run_dir(:)

  logical                                        :: times_medium
  logical                                        :: save_wfkern
  logical                                        :: save_snaps
  logical                                        :: norm_kernels
  logical                                        :: dump_vtk
  logical                                        :: dump_avs
  logical                                        :: dump_bin
  logical                                        :: dump_ascii
  logical                                        :: posteriori 

  integer, allocatable, dimension(:,:)           :: iidom_glob
  integer, allocatable, dimension(:)             :: idom_ipt
  integer                                        :: num_meshes, num_meshes_fwd
  integer, dimension(nproc_mesh)                 :: meshes, meshes_fwd
  integer, dimension(0:nproc_mesh-1)             :: countproc
  integer                                        :: ibwd_shift
  integer, allocatable, dimension(:)             :: npts_iproc 
  integer, allocatable, dimension(:)             :: npts_iproc_fwd  
  integer, allocatable, dimension(:)             :: n_src_iproc
  integer, allocatable, dimension(:,:)           :: iptkern_proc_rec
  integer, allocatable, dimension(:,:)           :: iptsem_proc_rec
  integer, allocatable, dimension(:,:)           :: iptkern_proc_src
  integer, allocatable, dimension(:,:)           :: iptsem_proc_src
  integer, allocatable, dimension(:,:)           :: iptkern_proc_srckern
  integer, allocatable, dimension(:)             :: iptsem_proc_srckern

  integer, allocatable, dimension(:,:)           :: iel_proc, iel_glob
  integer, allocatable, dimension(:,:,:,:,:)     :: gll_proc
  integer, allocatable, dimension(:,:)           :: imap1d

  integer, allocatable, dimension(:)             :: ind_cell2tot
  integer, allocatable, dimension(:)             :: valence_tot2cell
  integer, allocatable, dimension(:,:)           :: ind_tot2cell

  real(kind=realkind)                            :: phi0, dt, src_depth

  real(kind=realkind), dimension(ndumps,3)       :: u0, v0
  real(kind=realkind), allocatable               :: denomu(:), denomv(:)

  real(kind=realkind), allocatable, dimension(:) :: rho,lam,mu
  real(kind=realkind), allocatable, dimension(:) :: vp_fac, vs_fac 

  real(kind=realkind), dimension(ndisc)          :: discont
  real(kind=realkind)                            :: maxr, ar


  character(len=20)                              :: src_type1fwd
  character(len=20)                              :: src_type2fwd
  character(len=20)                              :: src_type1bwd
  character(len=20)                              :: src_type2bwd

  character(len=200)                             :: dir_fwdmesh
  character(len=200)                             :: dir_bwdmesh
  character(len=200)                             :: ext_mesh_name
  character(len=200), allocatable                :: dir_fwdmesh1(:)
  integer                                        :: lffwd, lfbwd, lfext

  logical                                        :: use_netcdf
  integer                                        :: ncid_fw_in, ncid_bw_in, &
                                                    ncid_fw_snapshot

  character(len=10),allocatable                  :: srctype(:)

  integer                                        :: ntw, nwin, num_it
  real(kind=realkind)                            :: twf1, twf2
  real(kind=realkind), allocatable, dimension(:) :: t1, t2
  integer, allocatable, dimension(:)             :: begdumps,enddumps
  integer, allocatable, dimension(:,:)           :: ind_win
  logical, allocatable, dimension(:,:)           :: inside_win, save_n_erase
  character(len=10), allocatable, dimension(:)   :: phase_name

  real(kind=realkind), allocatable               :: ext_mesh(:,:)
  integer                                        :: npts_tot, npts_cell
  real(kind=realkind), allocatable               :: xgd_tot(:), ygd_tot(:), zgd_tot(:)
  real(kind=realkind), allocatable               :: xgd_cell(:), ygd_cell(:), zgd_cell(:)
  real(kind=realkind), allocatable               :: xgd(:), ygd(:), zgd(:)
  real(kind=realkind), allocatable               :: phigd(:)
  real, allocatable                              :: W_vtk(:,:)
  real(kind=realkind), allocatable               :: phird(:)

  real(kind=realkind), allocatable               :: globmesh(:,:,:,:,:) 
                                            ! globmesh(ibeg:iend,ibeg:iend,nel,iproc,idim)
  
  real(kind=realkind), allocatable               :: mesh1d(:,:,:)
  real(kind=realkind), allocatable, dimension(:) :: azim1_fwd, azim2_fwd
  real(kind=realkind), allocatable, dimension(:) :: azim1_bwd, azim2_bwd
  real(kind=realkind), allocatable, dimension(:) :: azim1_src, azim2_src
  logical                                        :: do_azim,compute_src_kernels
  character(len=20)                              :: cell_mesh

  integer                                        :: npadd
  integer, dimension(:), allocatable             :: iproc_proc_tot
  integer, dimension(:), allocatable             :: iel_proc_tot,ipt_proc_tot

  real(kind=realkind)                            :: realjunk
  integer                                        :: intjunk
  character(len=50)                              :: charjunk,junk
  character(len=14)                              :: calc_or_load_kernels 
                                                    !calc_wfkernels or load_wfkernels

! SOURCE RELATED PARAMS
  integer                                        :: nsim
  real(kind=realkind)                            :: Mij(6),magnitude

! FROM SIMULATION.INFO OF SOLVER
  integer                                        :: ibeg_fwd, iend_fwd, ibeg_bwd, iend_bwd
  real(kind=realkind)                            :: shift_fact_fwd, shift_fact_deltat_fwd, &
                                                    shift_fact_seis_dt_fwd 
  real(kind=realkind)                            :: shift_fact_bwd, shift_fact_deltat_bwd, &
                                                    shift_fact_seis_dt_bwd 
  real(kind=realkind)                            :: shift_fact_deltat_coarse_fwd, &
                                                    shift_fact_deltat_coarse_bwd
  
! FOR CUBED SPHERE ONLY
  real(kind=realkind)                            :: ri   ! Inner spherical shell radius
  real(kind=realkind)                            :: re   ! Outer spherical shell radius
  integer                                        :: nang ! number of elements in the xi, 
                                                         ! eta direction
  integer                                        :: nr   ! number of element in the radial 
                                                         ! direction
  real(kind=realkind)                            :: dr, dang
  real(kind=realkind), dimension(:), allocatable :: xi_el, eta_el, r_el
  integer, dimension(:,:,:,:), allocatable       :: number_el
  real(kind=realkind), dimension(:,:,:,:), allocatable :: x_el, y_el, z_el
  integer :: nelt
  real(kind=realkind), dimension(:,:,:,:), allocatable :: xcol, ycol, zcol 
  real(kind=realkind), dimension(:,:,:,:), allocatable :: jacob   ! the jacobian
  real(kind=realkind), dimension(:,:,:,:), allocatable :: massmat ! the diagonal mass 
                                                                  ! matrix

  integer, dimension(:), allocatable             :: iel_cs, ipol_cs, jpol_cs, kpol_cs
  integer                                        :: nrpol, npol_cs
  real(kind=realkind), dimension(:), allocatable :: xi_i,xi_j,xi_k
  real(kind=realkind), dimension(:), allocatable :: wt_i,wt_j,wt_k

!=======================
end module data_mesh
!=======================
