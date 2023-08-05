#ifndef DP_FUNCTIONS_HPP
#define DP_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"
#include "physics/dp/dp_constants.hpp"

#include "share/scream_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream {
namespace dp {

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for DP. We use the ETI pattern for
 * these functions.
 *
 * DP assumptions:
 *  - Kokkos team policies have a vector length of 1
 */

struct element_t{};
struct hvcoord_t{};
struct timelevel_t{};
struct hybrid_t{};

template <typename ScalarT, typename DeviceT>
struct Functions
{
  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S>
  using BigPack = ekat::Pack<S,SCREAM_PACK_SIZE>;
  template <typename S>
  using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

  using IntSmallPack = SmallPack<Int>;
  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  using Mask  = ekat::Mask<Pack::n>;
  using Smask = ekat::Mask<Spack::n>;

  using KT = ekat::KokkosTypes<Device>;

  using C  = physics::Constants<Scalar>;
  using SC = dp::Constants<Scalar>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using view_3d = typename KT::template view_3d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

  template <typename S>
  using uview_2d = typename ekat::template Unmanaged<view_2d<S> >;

  using MemberType = typename KT::MemberType;

  using WorkspaceMgr = typename ekat::WorkspaceManager<Spack, Device>;
  using Workspace    = typename WorkspaceMgr::Workspace;

  //
  // --------- Functions ---------
  //

  KOKKOS_FUNCTION
  static void advance_iop_forcing(const Int& plev, const Int& pcnst, const Spack& scm_dt, const Spack& ps_in, const uview_1d<const Spack>& u_in, const uview_1d<const Spack>& v_in, const uview_1d<const Spack>& t_in, const uview_1d<const Spack>& q_in, const uview_1d<const Spack>& t_phys_frc, const uview_1d<Spack>& u_update, const uview_1d<Spack>& v_update, const uview_1d<Spack>& t_update, const uview_1d<Spack>& q_update);
  KOKKOS_FUNCTION
  static void advance_iop_nudging(const Int& plev, const Spack& scm_dt, const Spack& ps_in, const uview_1d<const Spack>& t_in, const uview_1d<const Spack>& q_in, const uview_1d<Spack>& t_update, const uview_1d<Spack>& q_update, const uview_1d<Spack>& relaxt, const uview_1d<Spack>& relaxq);
  KOKKOS_FUNCTION
  static void advance_iop_subsidence(const Int& plev, const Int& pcnst, const Spack& scm_dt, const Spack& ps_in, const uview_1d<const Spack>& u_in, const uview_1d<const Spack>& v_in, const uview_1d<const Spack>& t_in, const uview_1d<const Spack>& q_in, const uview_1d<Spack>& u_update, const uview_1d<Spack>& v_update, const uview_1d<Spack>& t_update, const uview_1d<Spack>& q_update);
  KOKKOS_FUNCTION
  static void iop_setinitial(const Int& nelemd, const uview_1d<element_t>& elem);
  KOKKOS_FUNCTION
  static void iop_broadcast();
  KOKKOS_FUNCTION
  static void apply_iop_forcing(const Int& nelemd, const uview_1d<element_t>& elem, hvcoord_t& hvcoord, const hybrid_t& hybrid, const timelevel_t& tl, const Int& n, const bool& t_before_advance, const Int& nets, const Int& nete);
  KOKKOS_FUNCTION
  static void iop_domain_relaxation(const Int& nelemd, const Int& np, const Int& nlev, const uview_1d<element_t>& elem, const hvcoord_t& hvcoord, const hybrid_t& hybrid, const Int& t1, const uview_1d<Spack>& dp, const Int& nelemd_todo, const Int& np_todo, const Spack& dt);
  KOKKOS_FUNCTION
  static void crm_resolved_turb(const Int& nelemd, const uview_1d<element_t>& elem, const hvcoord_t& hvcoord, const hybrid_t& hybrid, const Int& t1, const Int& nelemd_todo, const Int& np_todo);
  static void iop_default_opts(Spack& scmlat_out, Spack& scmlon_out, std::string& iopfile_out, bool& single_column_out, bool& scm_iop_srf_prop_out, bool& iop_nudge_tq_out, bool& iop_nudge_uv_out, Spack& iop_nudge_tq_low_out, Spack& iop_nudge_tq_high_out, Spack& iop_nudge_tscale_out, bool& scm_observed_aero_out, bool& iop_dosubsidence_out, bool& scm_multcols_out, bool& dp_crm_out, Spack& iop_perturb_high_out, bool& precip_off_out, bool& scm_zero_non_iop_tracers_out);
  static void iop_setopts(const Spack& scmlat_in, const Spack& scmlon_in, const std::string& iopfile_in, const bool& single_column_in, const bool& scm_iop_srf_prop_in, const bool& iop_nudge_tq_in, const bool& iop_nudge_uv_in, const Spack& iop_nudge_tq_low_in, const Spack& iop_nudge_tq_high_in, const Spack& iop_nudge_tscale_in, const bool& scm_observed_aero_in, const bool& iop_dosubsidence_in, const bool& scm_multcols_in, const bool& dp_crm_in, const Spack& iop_perturb_high_in, const bool& precip_off_in, const bool& scm_zero_non_iop_tracers_in);
  KOKKOS_FUNCTION
  static void setiopupdate_init();
  KOKKOS_FUNCTION
  static void setiopupdate();
  KOKKOS_FUNCTION
  static void readiopdata(const Int& plev, const bool& iop_update_phase1, const uview_1d<const Spack>& hyam, const uview_1d<const Spack>& hybm);
  KOKKOS_FUNCTION
  static void iop_intht();
}; // struct Functions

} // namespace dp
} // namespace scream

// If a GPU build, without relocatable device code enabled, make all code available
// to the translation unit; otherwise, ETI is used.
#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)  \
                                && !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)

# include "impl/dp_advance_iop_forcing_impl.hpp"
# include "impl/dp_advance_iop_nudging_impl.hpp"
# include "impl/dp_advance_iop_subsidence_impl.hpp"
# include "impl/dp_iop_setinitial_impl.hpp"
# include "impl/dp_iop_broadcast_impl.hpp"
# include "impl/dp_apply_iop_forcing_impl.hpp"
# include "impl/dp_iop_domain_relaxation_impl.hpp"
# include "impl/dp_crm_resolved_turb_impl.hpp"
# include "impl/dp_iop_default_opts_impl.hpp"
# include "impl/dp_iop_setopts_impl.hpp"
# include "impl/dp_setiopupdate_init_impl.hpp"
# include "impl/dp_setiopupdate_impl.hpp"
# include "impl/dp_readiopdata_impl.hpp"
# include "impl/dp_iop_intht_impl.hpp"
#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE

#endif // DP_FUNCTIONS_HPP
