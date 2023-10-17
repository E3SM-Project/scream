module hetfrz_classnuc_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for hetfrz_classnuc module.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver, begchunk, endchunk
use physconst,      only: rair, cpair, rh2o, rhoh2o, mwh2o, tmelt, pi
use constituents,   only: cnst_get_ind
use physics_types,  only: physics_state
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,   only: phys_getopts, use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_idx, rad_cnst_get_spec_idx, &
                            rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num, rad_cnst_get_mode_props

use physics_buffer, only: pbuf_add_field, dtype_r8, pbuf_old_tim_idx, &
                          pbuf_get_index, pbuf_get_field
use cam_history,    only: addfld, add_default, outfld

use ref_pres,       only: top_lev => trop_cloud_top_lev
use wv_saturation,  only: svp_water, svp_ice

use cam_logfile,    only: iulog
use error_messages, only: handle_errmsg, alloc_err
use cam_abortutils, only: endrun

use hetfrz_classnuc,   only: hetfrz_classnuc_init, hetfrz_classnuc_calc

implicit none
private
save

public :: &
   hetfrz_classnuc_cam_readnl,   &
   hetfrz_classnuc_cam_register
! Namelist variables
logical :: hist_hetfrz_classnuc = .false.

! Vars set via init method.
real(r8) :: mincld      ! minimum allowed cloud fraction

! constituent indices
integer :: &
   cldliq_idx = -1, &
   cldice_idx = -1, &
   numliq_idx = -1, &
   numice_idx = -1

! pbuf indices for fields provided by heterogeneous freezing
integer :: &
   frzimm_idx, &
   frzcnt_idx, &
   frzdep_idx

! pbuf indices for fields needed by heterogeneous freezing
integer :: &
   ast_idx = -1

! modal aerosols
integer, parameter :: MAM3_nmodes = 3
integer, parameter :: MAM7_nmodes = 7
integer, parameter :: MAM4_nmodes = 4
integer, parameter :: MAM5_nmodes = 5 
integer :: nmodes = -1             ! number of aerosol modes

! mode indices
integer :: mode_accum_idx  = -1    ! accumulation mode
integer :: mode_coarse_idx = -1    ! coarse mode
integer :: mode_finedust_idx = -1  ! fine dust mode
integer :: mode_coardust_idx = -1  ! coarse dust mode
integer :: mode_pcarbon_idx = -1   ! primary carbon mode

! mode properties
real(r8) :: alnsg_mode_accum
real(r8) :: alnsg_mode_coarse
real(r8) :: alnsg_mode_finedust
real(r8) :: alnsg_mode_coardust
real(r8) :: alnsg_mode_pcarbon

! specie properties
real(r8) :: specdens_dust
real(r8) :: specdens_so4
real(r8) :: specdens_bc
real(r8) :: specdens_soa
real(r8) :: specdens_pom
real(r8) :: specdens_mom

! List all species
integer :: ncnst = 0     ! Total number of constituents (mass and number) needed
                         ! by the parameterization (depends on aerosol model used)

integer :: so4_accum     ! sulfate in accumulation mode
integer :: bc_accum      ! black-c in accumulation mode
integer :: pom_accum     ! p-organic in accumulation mode
integer :: soa_accum     ! s-organic in accumulation mode
integer :: dst_accum     ! dust in accumulation mode
integer :: ncl_accum     ! seasalt in accumulation mode
integer :: mom_accum     ! marine-organic in accumulation mode
integer :: num_accum     ! number in accumulation mode

integer :: dst_coarse    ! dust in coarse mode
integer :: ncl_coarse    ! seasalt in coarse mode
integer :: so4_coarse    ! sulfate in coarse mode
integer :: bc_coarse     ! bc in coarse mode
integer :: pom_coarse    ! pom in coarse mode
integer :: soa_coarse    ! soa in coarse mode
integer :: mom_coarse    ! mom in coarse mode
integer :: num_coarse    ! number in coarse mode

integer :: dst_finedust  ! dust in finedust mode
integer :: so4_finedust  ! sulfate in finedust mode
integer :: num_finedust  ! number in finedust mode

integer :: dst_coardust  ! dust in coardust mode
integer :: so4_coardust  ! sulfate in coardust mode
integer :: num_coardust  ! number in coardust mode

integer :: bc_pcarbon    ! black-c in primary carbon mode
integer :: pom_pcarbon   ! p-organic in primary carbon mode
integer :: mom_pcarbon   ! marine-organic in primary carbon mode
integer :: num_pcarbon   ! number in primary carbon mode

! Index arrays for looping over all constituents
integer, allocatable :: mode_idx(:)
integer, allocatable :: spec_idx(:)

! Copy of cloud borne aerosols before modification by droplet nucleation
! The basis is converted from mass to volume.
real(r8), allocatable :: aer_cb(:,:,:,:)

! Copy of interstitial aerosols with basis converted from mass to volume.
real(r8), allocatable :: aer(:,:,:,:)

!===============================================================================
contains
!===============================================================================

subroutine hetfrz_classnuc_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'hetfrz_classnuc_cam_readnl'

  namelist /hetfrz_classnuc_nl/ hist_hetfrz_classnuc

  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'hetfrz_classnuc_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, hetfrz_classnuc_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

  end if

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(hist_hetfrz_classnuc, 1, mpilog, 0, mpicom)
#endif

end subroutine hetfrz_classnuc_cam_readnl

!================================================================================================

subroutine hetfrz_classnuc_cam_register()

   if (.not. use_hetfrz_classnuc) return

   ! pbuf fields provided by hetfrz_classnuc
   call pbuf_add_field('FRZIMM', 'physpkg', dtype_r8, (/pcols,pver/), frzimm_idx)
   call pbuf_add_field('FRZCNT', 'physpkg', dtype_r8, (/pcols,pver/), frzcnt_idx)
   call pbuf_add_field('FRZDEP', 'physpkg', dtype_r8, (/pcols,pver/), frzdep_idx)

end subroutine hetfrz_classnuc_cam_register

!================================================================================================



end module hetfrz_classnuc_cam
