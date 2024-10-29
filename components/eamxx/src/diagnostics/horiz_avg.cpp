#include "diagnostics/horiz_avg.hpp"

namespace scream {

HorizAvgDiag::HorizAvgDiag(const ekat::Comm &comm,
                           const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {
  const auto &fname = m_params.get<std::string>("field_name");
  m_diag_name       = fname + "_horiz_avg";
}

void HorizAvgDiag::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  const auto &fn = m_params.get<std::string>("field_name");
  const auto &gn = m_params.get<std::string>("grid_name");
  const auto g   = grids_manager->get_grid("Physics");
  m_area         = g->get_geometry_data("area").get_view<const Real *>();
  add_field<Required>(fn, gn);
}

void HorizAvgDiag::initialize_impl(const RunType /*run_type*/) {
  const auto &f = get_fields_in().front();
  using namespace ShortFieldTagsNames;
  const auto &fid    = f.get_header().get_identifier();
  const auto &layout = fid.get_layout();

  EKAT_REQUIRE_MSG(layout.rank() >= 1 && layout.rank() <= 4,
                   "Error! Field rank not supported by HorizAvgDiag.\n"
                   " - field name: " +
                       fid.name() +
                       "\n"
                       " - field layout: " +
                       layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.tags()[0] == COL,
                   "Error! HorizAvgDiag diagnostic expects a layout starting "
                   "with the 'COL' tag.\n"
                   " - field name  : " +
                       fid.name() +
                       "\n"
                       " - field layout: " +
                       layout.to_string() + "\n");

  FieldIdentifier d_fid(m_diag_name, layout.clone().strip_dim(COL),
                        fid.get_units(), fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  // get the area field
  int dim0     = layout.dim(0);
  const auto a = m_area;
  // calculate total area
  // m_total_area = 0.0;
  Kokkos::parallel_reduce(
      "HorizAvgDiag::compute_diagnostic_impl::total_area", dim0,
      KOKKOS_LAMBDA(const int icol, Real &accum) { accum += a(icol); }, m_total_area);
  // sum up the total area across ranks
  m_comm.all_reduce(&m_total_area, dim0, MPI_SUM);
  // ensure m_total_area is not zero
  m_total_area = m_total_area == 0.0 ? 1.0 : m_total_area;
}

void HorizAvgDiag::compute_diagnostic_impl() {
  using RangePolicy = Kokkos::RangePolicy<Field::device_t::execution_space>;
  using TeamPolicy  = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember  = typename TeamPolicy::member_type;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  const auto &f = get_fields_in().front();
  const auto &d = m_diagnostic_output;

  const auto &layout = f.get_header().get_identifier().get_layout();
  int dim0           = layout.dim(0);

  d.deep_copy(0);

  const auto a = m_area;
  const auto atot = m_total_area;

  switch(layout.rank()) {
    case 1: {
      auto f_view = f.get_view<const Real *>();
      auto d_view = d.get_view<Real>();

      auto p = ESU::get_default_team_policy(1, dim0);
      Kokkos::parallel_for(
          d.name(), p, KOKKOS_LAMBDA(const TeamMember &m) {
            Real sum = 0.0;
            Kokkos::parallel_reduce(
                Kokkos::TeamThreadRange(m, dim0),
                [&](const int icol, Real &accum) {
                  accum += (a(icol) / atot) * f_view(icol);
                },
                sum);
            Kokkos::single(Kokkos::PerTeam(m), [&]() { d_view() = sum; });
          });
    } break;
    case 2: {
      auto f_view = f.get_view<const Real **>();
      auto d_view = d.get_view<Real *>();

      const int dim1 = layout.dim(1);
      auto p         = ESU::get_default_team_policy(dim1, dim0);
      Kokkos::parallel_for(
          d.name(), p, KOKKOS_LAMBDA(const TeamMember &m) {
            const int j = m.league_rank();
            Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange(m, dim0),
                [&](int icol, Real &accum) {
                  accum += (a(icol) / atot) * f_view(icol, j);
                },
                d_view(j));
          });
    } break;
    case 3: {
      auto f_view = f.get_view<const Real ***>();
      auto d_view = d.get_view<Real **>();

      const int dim1 = layout.dim(1);
      const int dim2 = layout.dim(2);
      auto p         = ESU::get_default_team_policy(dim1 * dim2, dim0);
      Kokkos::parallel_for(
          d.name(), p, KOKKOS_LAMBDA(const TeamMember &m) {
            const int idx = m.league_rank();
            const int j   = idx / dim2;
            const int k   = idx % dim2;
            Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange(m, dim0),
                [&](int icol, Real &accum) {
                  accum += (a(icol) / atot) * f_view(icol, j, k);
                },
                d_view(j, k));
          });
    } break;
    case 4: {
      auto f_view = f.get_view<const Real ****>();
      auto d_view = d.get_view<Real ***>();

      const int dim1 = layout.dim(1);
      const int dim2 = layout.dim(2);
      const int dim3 = layout.dim(3);
      auto p         = ESU::get_default_team_policy(dim1 * dim2 * dim3, dim0);
      Kokkos::parallel_for(
          d.name(), p, KOKKOS_LAMBDA(const TeamMember &m) {
            const int idx = m.league_rank();
            const int j   = (idx / dim3) / dim2;
            const int k   = (idx / dim3) % dim2;
            const int l   = idx % dim3;
            Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange(m, dim0),
                [&](int icol, Real &accum) {
                  accum += (a(icol) / atot) * f_view(icol, j, k, l);
                },
                d_view(j, k, l));
          });
    } break;
  }
  Kokkos::fence();

#if SCREAM_MPI_ON_DEVICE
  m_comm.all_reduce(d.get_internal_view_data<Real>(), layout.size() / dim0,
                    MPI_SUM);
#else
  d.sync_to_host();
  m_comm.all_reduce(d.get_internal_view_data<Real, Host>(),
                    layout.size() / dim0, MPI_SUM);
  d.sync_to_dev();
#endif
}

}  // namespace scream
