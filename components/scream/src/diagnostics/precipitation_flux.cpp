#include "diagnostics/precipitation_rate.hpp"

namespace scream
{

// =========================================================================================
PrecipitationRateDiagnostic::RainWaterPathDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void PrecipitationRateDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto Q = kg/kg;
  Q.set_string("kg/kg");

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  FieldLayout scalar2d_layout_mid { {COL},     {m_num_cols}            };
  constexpr int ps = Pack::n;

  // The fields required for this diagnostic to be computed
  add_field<Required>("precip_liq_surf", scalar3d_layout_mid, Q, grid_name, ps);
  add_field<Required>("precip_ice_surf", scalar3d_layout_mid, Q, grid_name, ps);
  add_field<Required>("dt",             scalar3d_layout_mid, "s",  grid_name, ps);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar2d_layout_mid, m, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation();
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void PrecipitationRateDiagnostic::compute_diagnostic_impl()
{

  using PC         = scream::physics::Constants<Real>;
  constexpr Real gravit = PC::gravit;
  const auto npacks         = ekat::npack<Pack>(m_num_levs);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, npacks);
  const auto& precip_rate                = m_diagnostic_output.get_view<Real*>();
  const auto& precip_liq_surf             = get_field_in("precip_liq_surf").get_view<const Pack**>();
  const auto& precip_ice_surf             = get_field_in("precip_ice_surf").get_view<const Pack**>();\

  Kokkos::parallel_for("PrecipitationRateDiagnostic",
                       default_policy,
                       KOKKOS_LAMBDA(const MemberType& team) {
    const int icol = team.league_rank(); 
    precip_rate(icol) = (precip_liq_surf(icol) + precip_ice_surf(icol))/dt
  });
  Kokkos::fence();

  const auto ts = get_field_in("precip_liq_surf").get_header().get_tracking().get_time_stamp();
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);
}
// =========================================================================================
} //namespace scream
