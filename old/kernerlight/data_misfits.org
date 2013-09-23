  module data_misfits
  
  use global_parameters, only : realkind
  implicit none
  public 

  include 'mesh_params_kernel.h'

  logical :: consider_data
  character(len=100) :: data_dir
  real(kind=realkind), allocatable :: timewin(:),misfit_fct(:,:),chi(:,:)
  real(kind=realkind), allocatable :: usyn_filt(:,:),vsyn_filt(:,:),udat_filt(:,:)
  integer           :: num_filt,num_dir,num_misfits,num_tot_picks
  real(kind=realkind) :: time_data(1:ndumps),displ_data(1:ndumps)

  end module data_misfits
