module nucleate_ice_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for nucleate_ice module.
!
!  B. Eaton - Sept 2014
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver
use physconst,      only: pi, rair, tmelt
use constituents,   only: cnst_get_ind
use physics_types,  only: physics_state
use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,   only: use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num, rad_cnst_get_mode_props

use physics_buffer, only: pbuf_add_field, dtype_r8, pbuf_old_tim_idx, &
                          pbuf_get_index, pbuf_get_field
use cam_history,    only: addfld, add_default, outfld

use ref_pres,       only: top_lev => trop_cloud_top_lev
use wv_saturation,  only: qsat_water
#ifndef HAVE_ERF_INTRINSICS
use shr_spfn_mod,   only: erf => shr_spfn_erf
#endif
use cam_logfile,    only: iulog
use cam_abortutils, only: endrun

use nucleate_ice,   only: nucleati_init, nucleati


implicit none
private
save

public :: &
   nucleate_ice_cam_readnl,   &
   nucleate_ice_cam_register
   

! Namelist variables
logical, public, protected :: use_preexisting_ice = .false.
logical                    :: hist_preexisting_ice = .false.
logical, public, protected :: use_nie_nucleate = .false.
logical, public, protected :: use_dem_nucleate = .false.
real(r8)                   :: nucleate_ice_subgrid
real(r8)                   :: so4_sz_thresh_icenuc = huge(1.0_r8) !ice nucleation SO2 size threshold for aitken mode

! Vars set via init method.
real(r8) :: mincld      ! minimum allowed cloud fraction
real(r8) :: bulk_scale  ! prescribed aerosol bulk sulfur scale factor

! constituent indices
integer :: &
   cldliq_idx = -1, &
   cldice_idx = -1, &
   numice_idx = -1

integer :: &
   naai_idx,     &
   naai_hom_idx

integer :: &
   ast_idx   = -1, &
   dgnum_idx = -1

! Bulk aerosols
character(len=20), allocatable :: aername(:)
real(r8), allocatable :: num_to_mass_aer(:)

integer :: naer_all      ! number of aerosols affecting climate
integer :: idxsul   = -1 ! index in aerosol list for sulfate
integer :: idxdst1  = -1 ! index in aerosol list for dust1
integer :: idxdst2  = -1 ! index in aerosol list for dust2
integer :: idxdst3  = -1 ! index in aerosol list for dust3
integer :: idxdst4  = -1 ! index in aerosol list for dust4
integer :: idxbcphi = -1 ! index in aerosol list for Soot (BCPHIL)

! modal aerosols
logical :: clim_modal_aero
logical :: prog_modal_aero
real(r8) :: sigmag_accum
!logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
integer :: cnum_coarse_idx, ccoarse_dst_idx, ccoarse_so4_idx
integer :: cnum_strat_coarse_idx, cstrat_coarse_so4_idx
integer :: accum_so4_idx,accum_pom_idx,accum_soa_idx,accum_bc_idx,accum_dst_idx,accum_ncl_idx,accum_mom_idx

integer :: nmodes = -1
integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode

integer :: coarse_so4_idx = -1  ! index of so4 in coarse mode
! for so4 only aitken mode so4 to get so4_num for nucleation
! if defined MAM7S add stratosphere aitken to troposphere aitken
integer :: mode_strat_sulfate1_idx = -1
integer :: mode_strat_coarse_idx = -1
integer :: strat_coarse_so4_idx = -1
real(r8) :: sigmag_coarse, sigmag_strat_coarse
real(r8) :: sigmag_aitken
#if (defined MODAL_AERO_4MODE_MOM || defined MODAL_AERO_5MODE)
integer :: coarse_mom_idx = -1  ! index of mom in coarse mode
#endif

#if (defined RAIN_EVAP_TO_COARSE_AERO) 
integer :: coarse_bc_idx = -1  ! index of bc in coarse mode
integer :: coarse_pom_idx = -1  ! index of pom in coarse mode
integer :: coarse_soa_idx = -1  ! index of soa in coarse mode
#endif

integer :: mode_fine_dst_idx = -1   ! index of dust in fine dust mode
integer :: mode_pcarbon_idx  = -1  ! index of dust in accum mode
integer :: accum_dust_idx    = -1  ! index of dust in accum mode
integer :: fine_dust_idx    = -1   ! index of dust in fine mode

logical  :: separate_dust = .false.


!===============================================================================
contains
!===============================================================================

subroutine nucleate_ice_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'nucleate_ice_cam_readnl'

  namelist /nucleate_ice_nl/ use_preexisting_ice, hist_preexisting_ice, &
                             use_nie_nucleate, use_dem_nucleate,        &
                             nucleate_ice_subgrid, so4_sz_thresh_icenuc 

  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'nucleate_ice_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, nucleate_ice_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

  end if

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(use_preexisting_ice,  1, mpilog, 0, mpicom)
  call mpibcast(hist_preexisting_ice, 1, mpilog, 0, mpicom)
  call mpibcast(use_nie_nucleate, 1, mpilog, 0, mpicom)
  call mpibcast(use_dem_nucleate, 1, mpilog, 0, mpicom)
  call mpibcast(nucleate_ice_subgrid, 1, mpir8, 0, mpicom)
  call mpibcast(so4_sz_thresh_icenuc, 1, mpir8, 0, mpicom)
#endif

end subroutine nucleate_ice_cam_readnl

!================================================================================================

subroutine nucleate_ice_cam_register()

   call pbuf_add_field('NAAI',     'physpkg', dtype_r8, (/pcols,pver/), naai_idx)
   call pbuf_add_field('NAAI_HOM', 'physpkg', dtype_r8, (/pcols,pver/), naai_hom_idx)

end subroutine nucleate_ice_cam_register

end module nucleate_ice_cam
