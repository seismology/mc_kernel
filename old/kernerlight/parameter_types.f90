module parameter_types

  type kernelparamtype
      logical                           :: do_rho
      logical                           :: do_lam
      logical                           :: do_mu
      logical                           :: do_vp
      logical                           :: do_vs
      logical                           :: do_imped
      logical                           :: compute_src_kernels
      integer                           :: ntw
      real, allocatable                 :: t1(:), t2(:)
      character(len=10), allocatable    :: phase_name(:)
      character(len=4)                  :: char_wfkern
      real                              :: thr, phr, phi0
      real                              :: twf1, twf2
      real, allocatable                 :: tf1(:), tf2(:)
      logical                           :: times_medium, save_snaps, norm_kernels
      integer                           :: reccomp, dim_fwd
  end type

  type sourceparamtype
      character(len=20)                 :: type1
      character(len=20)                 :: type2
      character(len=10), allocatable    :: srctype(:)
      character(len=10), allocatable    :: exctype(:)
      character(len=20)                 :: stf_type
      integer                           :: nsim_req
      integer                           :: nsim
      real, dimension(6)                :: mij
      real, dimension(3)                :: forces
      real                              :: magnitude
      real                              :: depth
      real                              :: shift_fact
      real                              :: shift_fact_deltat
      real                              :: shift_fact_seis_dt
      real                              :: shift_fact_deltat_coarse
  end type

  type filterparamtype
      complex(kind=8),allocatable       :: specfilt(:)
      real(kind=8)                      :: power
      integer                           :: what
      character(len=3)                  :: type
      real(kind=8)                      :: period, period_low, period_hi
      real(kind=8)                      :: fmin, fmax
      real(kind=8)                      :: omegamin, omegamax
  end type

  type ncparamtype
      integer                           :: ncid
      integer                           :: snap, surf, mesh
      character(len=200)                :: meshdir
      integer                           :: ndumps
      logical                           :: ordered_output
  end type

  type ioparamtype
      character(len=200)                :: dir_fwdmesh
      character(len=200)                :: dir_bwdmesh
      character(len=200)                :: ext_mesh_name
      character(len=200)                :: data_dir
      character(len=200), allocatable   :: dir_fwdmesh1(:)
      integer                           :: lffwd, lfbwd, lfext
      integer                           :: nsim_fwd, nsim_bwd
      type(ncparamtype), allocatable    :: fwd(:)
      type(ncparamtype), allocatable    :: bwd(:)
  end type

  type paramtype
      type(kernelparamtype) :: k   !< Kernel
      type(sourceparamtype) :: src_fwd, src_bwd !< Source
      type(ioparamtype)     :: io  !< I/O
      type(filterparamtype) :: filter !< Filter
  end type
    

end module parameter_types


