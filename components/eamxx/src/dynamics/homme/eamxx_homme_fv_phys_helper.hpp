#ifndef EAMXX_HOMME_FV_PHYS_HELPER_HPP
#define EAMXX_HOMME_FV_PHYS_HELPER_HPP

#include "dynamics/homme/homme_dimensions.hpp"

#include "share/grid/grids_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/field/field.hpp"

#include <ekat/kokkos/ekat_kokkos_types.hpp>
#include <ekat/ekat_pack.hpp>

namespace scream {

struct HommeFvPhysHelper {
public:
  using KT = KokkosTypes<DefaultDevice>;
  template<typename T,int N>
  using view_ND = typename KT::view_ND<T,N>;

  static constexpr int N = HOMMEXX_PACK_SIZE;
  using PackT = ekat::Pack<Real,N>;

  static HommeFvPhysHelper& instance () {
    static HommeFvPhysHelper dph;
    return dph;
  }

  void set_grids (const std::shared_ptr<const AbstractGrid>& dyn_grid,
                  const std::shared_ptr<const AbstractGrid>& cgll_grid,
                  const std::shared_ptr<const AbstractGrid>& phys_grid);

  void requested_buffer_size_in_bytes ();

  void dyn_to_fv_phys_init (const bool restart, const util::TimeStamp& t0);

  void remap_dyn_to_fv_phys () const;
  void remap_dyn_to_fv_phys (const view_ND<Real,2>& T_out, const view_ND<Real,3>& uv_out, const view_ND<Real,3>& Q_out) const;
  void remap_fv_phys_to_dyn () const;

  void clean_up ();

  bool fv_phys_active = false;
  int  pgN            = -1;

  // Copy physics T,uv state to FT,M to form tendencies in next dynamics step.
  static void copy_prev (const int ncols, const int nlevs,
                         const view_ND<PackT,2>& T,  const view_ND<PackT,3>& uv,
                         const view_ND<PackT,2>& FT, const view_ND<PackT,3>& FM);

  // Physics forcing
  Field m_FT_phys;
  Field m_FM_phys;

  // Dynamics state on physics grid
  Field m_T_phys;
  Field m_uv_phys;
  Field m_ps_phys;
  Field m_phis_phys;
  Field m_omega_phys;
  Field m_dp_phys;
  Field m_Q_phys;

  // During restart, we can't overwrite T_mid, horiz_winds, and tracers, as these
  // contain data used in homme_pre_process to compute tendencies.

  std::shared_ptr<const AbstractGrid> m_dyn_grid;
  std::shared_ptr<const AbstractGrid> m_cgll_grid;
  std::shared_ptr<const AbstractGrid> m_phys_grid;

  HommeFvPhysHelper () = default;
};

} // namespace scream

#endif  // EAMXX_HOMME_FV_PHYS_HELPER_HPP
