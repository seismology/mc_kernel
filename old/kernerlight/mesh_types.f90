module mesh_types
    type meshtype
        integer, allocatable, dimension(:)   :: idom
        real, allocatable, dimension(:)      :: s, z
        integer                              :: npoints
    end type

    type modeltype
        real, allocatable, dimension(:)      :: discs
        integer                              :: ndisc
        real                                 :: maxr
        character(len=200)                   :: background_model
    end type

    type kernelmeshtype
        integer, allocatable, dimension(:)   :: idom
        integer, allocatable, dimension(:)   :: idom_all
        real, allocatable, dimension(:)      :: x, y, z
        real, allocatable, dimension(:)      :: x_all, y_all, z_all
        real, allocatable, dimension(:)      :: fwd_pointid
        real, allocatable, dimension(:)      :: bwd_pointid
 !       real, dimension(:), pointer                :: x_now, y_now, z_now
 !       real, dimension(:), pointer                :: fwd_id_now, bwd_id_now
        integer                              :: npoints, npoints_all
    end type

    type kerneltype
        real(kind=4), allocatable, dimension(:) :: lamkernel
        real(kind=4), allocatable, dimension(:) :: lamttkern,lamampkern
        real(kind=4), allocatable, dimension(:) :: mukernel
        real(kind=4), allocatable, dimension(:) :: muttkern,muampkern
        real(kind=4), allocatable, dimension(:) :: rhokernel
        real(kind=4), allocatable, dimension(:) :: rhottkern,rhoampkern
        real(kind=4), allocatable, dimension(:) :: vpkernel
        real(kind=4), allocatable, dimension(:) :: vpttkern,vpampkern
        real(kind=4), allocatable, dimension(:) :: vskernel
        real(kind=4), allocatable, dimension(:) :: vsttkern,vsampkern
        real(kind=4), allocatable, dimension(:) :: impkernel
        real(kind=4), allocatable, dimension(:) :: impttkern,impampkern

! The expensive folks: space-time/frequency arrays
        real(kind=4), allocatable, dimension(:,:) :: lamkern_kxt,mukern_kxt,rhokern_kxt
        real(kind=4), allocatable, dimension(:,:,:) :: field_phys_dim
      
        real(kind=8), allocatable, dimension(:,:) :: kerphys
        complex(kind=8), allocatable, dimension(:,:,:) :: ugd_src_spec,ugd_rec_spec
        complex(kind=8), allocatable, dimension(:,:) :: kerspec
    end type

end module
