#include "diagnostics/meridional_wind.hpp"

namespace scream
{

// =========================================================================================
MeridionalWindDiagnostic::MeridionalWindDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void MeridionalWindDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto vel = m/s;
  vel.set_string("m/s");

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  FieldLayout scalar2d_layout_mid { {COL},     {m_num_cols}            };
  FieldLayout horiz_wind_layout { {COL,CMP,LEV}, {m_num_cols,2,m_num_levs} };
  constexpr int ps = Pack::n;

  // The fields required for this diagnostic to be computed
  // Note both u and v are packaged into the single field horiz_winds.
  add_field<Required>("horiz_winds",    horiz_wind_layout,   m/s, grid_name, ps);


  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar2d_layout_mid, m, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation();
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void MeridionalWindDiagnostic::compute_diagnostic_impl()
{
  const auto npacks  = ekat::npack<Pack>(m_num_levs);
  const auto& horiz_winds        = get_field_in("horiz_winds").get_view<const Pack***>();
  const auto& v = m_diagnostic_output.get_view<Pack**>();

  Kokkos::parallel_for("MeridionalWindDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*npacks),
                       KOKKOS_LAMBDA (const int& idx) {
      const int icol  = idx / npacks;
      const int jpack = idx % npacks;
      v(icol,jpack) = horiz_winds(icol,1,jpack);
  });
  Kokkos::fence();


}
// =========================================================================================
} //namespace scream
