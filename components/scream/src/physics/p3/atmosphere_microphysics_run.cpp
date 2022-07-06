#include "physics/p3/atmosphere_microphysics.hpp"

namespace scream {

void P3Microphysics::run_impl (const int dt_proc)
{

  // TODO: Need to check if this is the first timestep of a subcycling, this will be taken care of
  // in a separate PR and then these lines can be updated.
  // We also need to have the dt_atm be explicitly passed to P3 somehow.  Will be taken care of in a separate PR.
  const auto& precip_liq_surf = get_field_out("precip_liq_surf").get_view<Real*>();
  const auto& precip_ice_surf = get_field_out("precip_ice_surf").get_view<Real*>();
  Real dt_atm = m_params.get<Real>("dt_atm",1800);
  Int  num_subcycles = m_params.get<Int>("dt_subcycle",6);
  if (m_subcycle_step%num_subcycles==0) {
    Kokkos::deep_copy(precip_liq_surf,0.0);
    Kokkos::deep_copy(precip_ice_surf,0.0);
    m_subcycle_step = 0;
  }
  m_subcycle_step += 1;
  // Assign values to local arrays used by P3, these are now stored in p3_loc.
  Kokkos::parallel_for(
    "p3_main_local_vals",
    Kokkos::RangePolicy<>(0,m_num_cols),
    p3_preproc
  ); // Kokkos::parallel_for(p3_main_local_vals)
  Kokkos::fence();

  // Update the variables in the p3 input structures with local values.

  infrastructure.dt = dt_proc;
  infrastructure.it++;

  // Reset internal WSM variables.
  workspace_mgr.reset_internals();

  // Run p3 main
  P3F::p3_main(prog_state, diag_inputs, diag_outputs, infrastructure,
               history_only, lookup_tables, workspace_mgr, m_num_cols, m_num_levs);

  // Rescale precip fluxes from P3 for the atmosphere timestep in case
  // there is subcycling.
  Kokkos::parallel_for("",m_num_cols, [&] (const int& icol) {
    diag_outputs.precip_liq_surf(icol) *= dt_proc/dt_atm;
    diag_outputs.precip_ice_surf(icol) *= dt_proc/dt_atm;
  });

  // Conduct the post-processing of the p3_main output.
  Kokkos::parallel_for(
    "p3_main_local_vals",
    Kokkos::RangePolicy<>(0,m_num_cols),
    p3_postproc
  ); // Kokkos::parallel_for(p3_main_local_vals)
  Kokkos::fence();
}

} // namespace scream
