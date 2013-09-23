  module data_rota
  
  use global_parameters, only : realkind
  implicit none

  public 
  integer, dimension(:), allocatable :: elem_d_src
  integer, dimension(:), allocatable :: iproc_d_src
!
  integer, dimension(:), allocatable :: elem_d_rec
  integer, dimension(:), allocatable :: iproc_d_rec

  integer, dimension(:), allocatable :: ipt_d_rec
  integer, dimension(:), allocatable :: ipt_d_src

  real(kind=realkind), dimension(:,:,:),allocatable :: RT,transRT,tens6x6_RT
  real(kind=realkind), dimension(3,3) :: vec_rot,trans_vec_rot
  real(kind=realkind), dimension(6,6) :: tens_rot
  real(kind=realkind), dimension(:,:), allocatable :: rot_cyl2cart
  real(kind=realkind), dimension(:),allocatable :: cosph,sinph
! source kernel
  integer, dimension(:), allocatable :: ipt_srckern,iproc_srckern,iel_srckern
  real(kind=realkind), dimension(:,:), allocatable :: s_srckern,z_srckern

  end module data_rota
