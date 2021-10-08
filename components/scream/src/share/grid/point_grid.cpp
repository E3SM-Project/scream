#include "share/grid/point_grid.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/share/physics_constants.hpp"
#include "ekat/mpi/ekat_comm.hpp"

#include <numeric>

namespace scream {

PointGrid::
PointGrid (const std::string& grid_name,
           const int          num_global_cols,
           const int          num_my_cols,
           const int          num_vertical_levels)
 : AbstractGrid(GridType::Point,grid_name)
{
  m_num_global_dofs = num_global_cols;
  m_num_local_dofs  = num_my_cols;
  m_num_vert_levs   = num_vertical_levels;

  // Sanity checks
  EKAT_REQUIRE_MSG (m_num_local_dofs>=0,
                    "Error! The number of local columns is negative.\n");
  EKAT_REQUIRE_MSG (num_vertical_levels>=1,
                    "Error! The number of vertical levels is not positive.\n");
  EKAT_REQUIRE_MSG (m_num_global_dofs>=m_num_local_dofs,
                    "Error! The number of global columns is smaller than the local one.\n");

  // The lid->idx map is the identity map.
  m_lid_to_idx = decltype(m_lid_to_idx)("lid to idx",m_num_local_dofs,1);

  auto h_lid_to_idx = Kokkos::create_mirror_view(m_lid_to_idx);
  std::iota(h_lid_to_idx.data(),h_lid_to_idx.data()+m_num_local_dofs,0);
  Kokkos::deep_copy(m_lid_to_idx,h_lid_to_idx);
}

FieldLayout
PointGrid::get_2d_scalar_layout () const
{
  using namespace ShortFieldTagsNames;

  return FieldLayout({COL},{m_num_local_dofs});
}

FieldLayout
PointGrid::get_2d_vector_layout (const FieldTag vector_tag, const int vector_dim) const
{
  using namespace ShortFieldTagsNames;

  return FieldLayout({COL,vector_tag},{m_num_local_dofs,vector_dim});
}

FieldLayout
PointGrid::get_3d_scalar_layout (const bool midpoints) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  return FieldLayout({COL,VL},{m_num_local_dofs,nvl});
}

FieldLayout
PointGrid::get_3d_vector_layout (const bool midpoints, const FieldTag vector_tag, const int vector_dim) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  return FieldLayout({COL,vector_tag,VL},{m_num_local_dofs,vector_dim,nvl});
}

void PointGrid::
set_dofs (const dofs_list_type& dofs)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (
    dofs.extent_int(0)==m_num_local_dofs,
    "Error! Input views have the wrong size:\n"
    "   expected: " + std::to_string(m_num_local_dofs) + "\n"
    "   provided: " + std::to_string(dofs.extent_int(0)) + "\n");

  m_dofs_gids  = dofs;
}

void PointGrid::
set_geometry_data (const std::string& name, const geo_view_type& data) {
  // Sanity checks
  EKAT_REQUIRE_MSG (data.extent_int(0)==m_num_local_dofs,
                    "Error! Input geometry data has wrong dimensions.\n");
  EKAT_REQUIRE_MSG (name=="lat" || name=="lon" || name=="area",
                    "Error! Point grid does not support geometry data '" + name + "'.\n");

  m_geo_views[name] = data;
}

std::shared_ptr<const PointGrid>
create_point_grid (const std::string& grid_name,
                   const int num_global_cols,
                   const int num_vertical_lev,
                   const ekat::Comm& comm)
{
  // Compute how many columns are owned by this rank
  const int num_procs = comm.size();

  auto num_my_cols = num_global_cols / num_procs;
  int remainder   = num_global_cols % num_procs;
  int dof_offset  = num_my_cols*comm.rank();
  if (comm.rank() < remainder) {
    ++num_my_cols;
    dof_offset += comm.rank();
  } else {
    dof_offset += remainder;
  }

  auto grid = std::make_shared<PointGrid>(grid_name,num_global_cols,num_my_cols,num_vertical_lev);

  PointGrid::dofs_list_type dofs_gids ("phys dofs",num_my_cols);
  auto h_dofs_gids = Kokkos::create_mirror_view(dofs_gids);
  std::iota(h_dofs_gids.data(),h_dofs_gids.data()+num_my_cols,dof_offset);
  Kokkos::deep_copy(dofs_gids,h_dofs_gids);

  grid->set_dofs(dofs_gids);

  using device_type       = DefaultDevice;
  using kokkos_types      = KokkosTypes<device_type>;
  using geo_view_type     = kokkos_types::view_1d<Real>;
  using KT                = KokkosTypes<DefaultDevice>;
  using C                 = scream::physics::Constants<Real>;

  // Store cell area, longitude, and latitude in geometry data.
  // For  longitude and latitude, set values to NaN since they
  // are currently not required from any application using PointGrid
  geo_view_type area("area", num_my_cols);
  geo_view_type lon ("lon",  num_my_cols);
  geo_view_type lat ("lat",  num_my_cols);

  // Estimate cell area for a uniform grid by taking the surface area
  // of the earth divided by the number of columns
  const Real rearth    = C::r_earth;
  const Real pi        = C::Pi;
  const Real cell_area = 4*pi*rearth*rearth/num_my_cols;

  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(num_my_cols, num_vertical_lev);
  Kokkos::parallel_for("area_loop", policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();
    area(i) = cell_area;
    lon(i) = std::nan("");
    lat(i) = std::nan("");
  });
  grid->set_geometry_data("area", area);
  grid->set_geometry_data("lon",  lon);
  grid->set_geometry_data("lat",  lat);

  return grid;
}

} // namespace scream
