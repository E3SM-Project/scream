#include "eamxx_homme_process_interface.hpp"

// HOMMEXX Includes
#include "Context.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"
#include "GllFvRemap.hpp"

// Scream includes
#include "share/util/eamxx_fv_phys_rrtmgp_active_gases_workaround.hpp"
#include "dynamics/homme/homme_dimensions.hpp"

namespace scream {

// See the [rrtmgp active gases] note in share/util/eamxx_fv_phys_rrtmgp_active_gases_workaround.hpp
void HommeDynamics
::fv_phys_rrtmgp_active_gases_init (const std::shared_ptr<const GridsManager>& gm) {
  auto& trace_gases_workaround = TraceGasesWorkaround::singleton();
  if (trace_gases_workaround.is_restart()) return; // always false b/c it hasn't been set yet
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  auto molmol = mol/mol;
  molmol.set_string("mol/mol");
  const auto& rgn = m_cgll_grid->name();
  const auto& pgn = m_phys_grid->name();
  const auto rnc = m_cgll_grid->get_num_local_dofs();
  const auto pnc = m_phys_grid->get_num_local_dofs();
  const auto nlev = m_cgll_grid->get_num_vertical_levels();
  for (const auto& e : trace_gases_workaround.get_active_gases()) {
    add_field<Required>(e, FieldLayout({COL,LEV},{rnc,nlev}), molmol, rgn);
    // 'Updated' rather than just 'Computed' so that it gets written to the
    // restart file.
    add_field<Updated >(e, FieldLayout({COL,LEV},{pnc,nlev}), molmol, pgn);
  }
  trace_gases_workaround.set_remapper(gm->create_remapper(m_cgll_grid, m_dyn_grid));
}

// See the [rrtmgp active gases] note in share/util/eamxx_fv_phys_rrtmgp_active_gases_workaround.hpp
void HommeDynamics::fv_phys_rrtmgp_active_gases_remap (const int pgN) {
  // Note re: restart: Ideally, we'd know if we're restarting before having to
  // call add_field above. However, we only find out after. Because the pg2
  // field was declared Updated, it will read the restart data. But we don't
  // actually want to remap from CGLL to pg2 now. So if restarting, just do the
  // cleanup part at the end.

  auto& trace_gases_workaround = TraceGasesWorkaround::singleton();
  const auto& rgn = m_cgll_grid->name();
  if (not trace_gases_workaround.is_restart()) {
    using namespace ShortFieldTagsNames;
    const auto& dgn = m_dyn_grid ->name();
    const auto& pgn = m_phys_grid->name();
    constexpr int NGP = HOMMEXX_NP;
    const int ngll = NGP*NGP;
    const int npg = pgN*pgN;
    const int nelem = m_dyn_grid->get_num_local_dofs()/ngll;
    { // CGLL -> DGLL
      const auto nlev = m_dyn_grid->get_num_vertical_levels();
      for (const auto& e : trace_gases_workaround.get_active_gases())
        create_helper_field(e, {EL,GP,GP,LEV}, {nelem,NGP,NGP,nlev}, dgn);
      auto r = trace_gases_workaround.get_remapper();
      r->registration_begins();
      for (const auto& e : trace_gases_workaround.get_active_gases())
        r->register_field(get_field_in(e, rgn), m_helper_fields.at(e));
      r->registration_ends();
      r->remap(true);
      trace_gases_workaround.erase_remapper();
    }
    { // DGLL -> PGN
      const auto& c = Homme::Context::singleton();
      auto& gfr = c.get<Homme::GllFvRemap>();
      const auto time_idx = c.get<Homme::TimeLevel>().n0;
      for (const auto& e : trace_gases_workaround.get_active_gases()) {
        const auto& f_dgll = m_helper_fields.at(e);
        const auto& f_phys = get_field_out(e, pgn);
        const auto& v_dgll = f_dgll.get_view<const Real****>();
        const auto& v_phys = f_phys.get_view<Real**>();
        assert(v_dgll.extent_int(0) == nelem and
               v_dgll.extent_int(1)*v_dgll.extent_int(2) == ngll);
        const auto in_dgll = Homme::GllFvRemap::CPhys3T(
          v_dgll.data(), nelem, 1, ngll, v_dgll.extent_int(3));
        assert(nelem*npg == v_phys.extent_int(0));
        const auto out_phys = Homme::GllFvRemap::Phys3T(
          v_phys.data(), nelem, npg, 1, v_phys.extent_int(1));
        gfr.remap_tracer_dyn_to_fv_phys(time_idx, 1, in_dgll, out_phys);
        Kokkos::fence();
      }
    }
  }
  // Done with all of these, so remove them.
  trace_gases_workaround.erase_remapper();
  for (const auto& e : trace_gases_workaround.get_active_gases())
    m_helper_fields.erase(e);
  for (const auto& e : trace_gases_workaround.get_active_gases())
    remove_field(e, rgn);
}

} // namespace scream
