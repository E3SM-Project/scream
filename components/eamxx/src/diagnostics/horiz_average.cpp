#include "diagnostics/horiz_average.hpp"

namespace scream
{

// =========================================================================================
HorizAverage::HorizAverage (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  const auto& fname = m_params.get<std::string>("field_name");
  m_diag_name = fname + "_havg";
}

void HorizAverage::
set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  const auto& fname = m_params.get<std::string>("field_name");
  const auto& gname = m_params.get<std::string>("grid_name");
  add_field<Required>(fname,gname);
}

void HorizAverage::
initialize_impl (const RunType /*run_type*/)
{
  const auto& f   = get_fields_in().front();
  // Sanity checks
  using namespace ShortFieldTagsNames;
  const auto& fid    = f.get_header().get_identifier();
  const auto& layout = fid.get_layout();

  EKAT_REQUIRE_MSG (layout.rank()>1 && layout.rank()<=4,
      "Error! Field rank not supported by HorizAverage.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");
  EKAT_REQUIRE_MSG (layout.tags()[0]==COL,
      "Error! HorizAverage diagnostic expects a layout starting with the 'COL' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // All good, create the diag output
  FieldIdentifier d_fid (m_diag_name,layout.clone().strip_dim(COL),fid.get_units(),fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();
}

// =========================================================================================
void HorizAverage::compute_diagnostic_impl()
{
  using RangePolicy = Kokkos::RangePolicy<Field::device_t::execution_space>;
  using TeamPolicy = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember = typename TeamPolicy::member_type;

  const auto& f = get_fields_in().front();
  const auto& d = m_diagnostic_output;
  d.deep_copy(0);

  const auto& layout = f.get_header().get_identifier().get_layout();
  int dim0 = layout.dim(0);
  switch (layout.rank()) {
    case 1:
      {
        // 1-dim is different than 2+.
        auto f_view = f.get_view<const Real*>();

        auto& result = d.get_view<Real,Host>()();
        RangePolicy p(0,dim0);
        Kokkos::parallel_reduce(d.name(),p,KOKKOS_LAMBDA(const int icol,Real& accum) {
          accum += f_view(icol);
        },d.get_view<Real,Host>()());
      }
    case 2:
      {
        auto f_view = f.get_view<const Real**>();
        auto d_view = d.get_view<      Real*>();

        const int dim1 = diag_layout.dim(1);
        TeamPolicy p(dim1,dim0);
        Kokkos::parallel_for(d.name(),p,KOKKOS_LAMBDA(const TeamMember& m) {
            const int j = m.league_rank();
            Kokkos::parallel_reduce(Kokkos::TeamVectorRange(0,dim0),[&](int icol,Real& accim){
              accum += f_view(icol,j);
            },d_view(j));
        });
      }
      break;
    case 3:
      {
        auto f_view = f.get_view<const Real***>();
        auto d_view = d.get_view<      Real**>();

        const int dim1 = diag_layout.dim(1);
        const int dim2 = diag_layout.dim(2);
        TeamPolicy p(dim1*dim2,dim0);
        Kokkos::parallel_for(d.name(),p,KOKKOS_LAMBDA(const TeamMember& m) {
            const int idx = m.league_rank();
            const int j = idx / dim2;
            const int k = idx % dim2;
            Kokkos::parallel_reduce(Kokkos::TeamVectorRange(0,dim0),[&](int icol,Real& accim){
              accum += f_view(icol,j,k);
            },d_view(j,k));
        });
      }
      break;
    case 4:
      {
        auto f_view = f.get_view<const Real****>();
        auto d_view = d.get_view<      Real***>();

        const int dim1 = diag_layout.dim(1);
        const int dim2 = diag_layout.dim(2);
        const int dim3 = diag_layout.dim(3);
        TeamPolicy p(dim1*dim2*dim3,dim0);
        Kokkos::parallel_for(d.name(),p,KOKKOS_LAMBDA(const TeamMember& m) {
            const int idx = m.league_rank();
            const int j = (idx / dim3) / dim2;
            const int k = (idx / dim3) % dim2;
            const int l =  idx % dim3;
            Kokkos::parallel_reduce(Kokkos::TeamVectorRange(0,dim0),[&](int icol,Real& accim){
              accum += f_view(icol,j,k);
            },d_view(j,k));
        });
      }
      break;
  }
  Kokkos::fence();

#if SCREAM_MPI_ON_DEVICE
  m_comm.all_reduce(d.get_internal_view_data<Real>,layout.size()/dim0,MPI_SUM);
#else
  d.sync_to_host();
  m_comm.all_reduce(d.get_internal_view_data<Real,Host>,layout.size()/dim0,MPI_SUM);
  d.sync_to_dev();
#endif
}

} // namespace scream
