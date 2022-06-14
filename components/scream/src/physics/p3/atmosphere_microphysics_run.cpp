#include "physics/p3/atmosphere_microphysics.hpp"

namespace scream {

void P3Microphysics::run_impl (const int dt)
{
  // Grab the fields that will record the water mass after p3
  const auto& qv_before_p3  = get_field_out("qv_before_p3").get_view<Pack**>();
  const auto& qc_before_p3  = get_field_out("qc_before_p3").get_view<Pack**>();
  const auto& qr_before_p3  = get_field_out("qr_before_p3").get_view<Pack**>();
  const auto& qi_before_p3  = get_field_out("qi_before_p3").get_view<Pack**>();

  Kokkos::deep_copy(qv_before_p3,prog_state.qv);
  Kokkos::deep_copy(qc_before_p3,prog_state.qc;
  Kokkos::deep_copy(qr_before_p3,prog_state.qr);
  Kokkos::deep_copy(qi_before_p3,prog_state.qi);
  // Assign values to local arrays used by P3, these are now stored in p3_loc.
  Kokkos::parallel_for(
    "p3_main_local_vals",
    Kokkos::RangePolicy<>(0,m_num_cols),
    p3_preproc
  ); // Kokkos::parallel_for(p3_main_local_vals)
  Kokkos::fence();

  // Update the variables in the p3 input structures with local values.

  infrastructure.dt = dt;
  infrastructure.it++;

  // Reset internal WSM variables.
  workspace_mgr.reset_internals();

  // Run p3 main
  P3F::p3_main(prog_state, diag_inputs, diag_outputs, infrastructure,
               history_only, lookup_tables, workspace_mgr, m_num_cols, m_num_levs);

  // Conduct the post-processing of the p3_main output.
  Kokkos::parallel_for(
    "p3_main_local_vals",
    Kokkos::RangePolicy<>(0,m_num_cols),
    p3_postproc
  ); // Kokkos::parallel_for(p3_main_local_vals)
  Kokkos::fence();

  // Grab the fields that will record the water mass after p3
  const auto& qv_after_p3  = get_field_out("qv_after_p3").get_view<Pack**>();
  const auto& qc_after_p3  = get_field_out("qc_after_p3").get_view<Pack**>();
  const auto& qr_after_p3  = get_field_out("qr_after_p3").get_view<Pack**>();
  const auto& qi_after_p3  = get_field_out("qi_after_p3").get_view<Pack**>();

  Kokkos::deep_copy(qv_after_p3,prog_state.qv);
  Kokkos::deep_copy(qc_after_p3,prog_state.qc;
  Kokkos::deep_copy(qr_after_p3,prog_state.qr);
  Kokkos::deep_copy(qi_after_p3,prog_state.qi);
}

} // namespace scream
