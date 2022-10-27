#include "share/grid/mesh_free_grids_manager.hpp"

#include "share/io/scorpio_input.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/remap/do_nothing_remapper.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

#include <memory>
#include <numeric>

namespace scream {

MeshFreeGridsManager::
MeshFreeGridsManager (const ekat::Comm& comm, const ekat::ParameterList& p)
 : m_params (p)
 , m_comm   (comm)
{
}

MeshFreeGridsManager::remapper_ptr_type
MeshFreeGridsManager::
do_create_remapper (const grid_ptr_type from_grid,
                    const grid_ptr_type to_grid) const
{
  return std::make_shared<DoNothingRemapper>(from_grid,to_grid);
}

void MeshFreeGridsManager::
build_grids ()
{
  auto has_positive_int = [&](const std::string& n) -> bool {
    return m_params.isParameter(n) && (m_params.get<int>(n)>0);
  };
  const bool build_pt = has_positive_int("number_of_global_columns");
  const bool build_se = has_positive_int("number_of_local_elements") &&
                        has_positive_int("number_of_gauss_points");

  const int num_vertical_levels = m_params.get<int>("number_of_vertical_levels");

  if (build_se) {
    // Build a set of completely disconnected spectral elements.
    const int num_local_elems  = m_params.get<int>("number_of_local_elements");
    const int num_gp           = m_params.get<int>("number_of_gauss_points");

    // Create the grid
    std::shared_ptr<SEGrid> se_grid;
    se_grid = std::make_shared<SEGrid>("SE Grid",num_local_elems,num_gp,num_vertical_levels,m_comm);
    se_grid->setSelfPointer(se_grid);

    // Set up the degrees of freedom.
    auto dof_gids = se_grid->get_dofs_gids();
    auto lid2idx  = se_grid->get_lid_to_idx_map();

    auto host_dofs    = dof_gids.template get_view<AbstractGrid::gid_type*,Host>();
    auto host_lid2idx = lid2idx.template get_view<int**,Host>();

    // Count unique local dofs. On all elems except the very last one (on rank N),
    // we have num_gp*(num_gp-1) unique dofs;
    int num_local_dofs = num_local_elems*num_gp*num_gp;
    int offset = num_local_dofs*m_comm.rank();

    for (int ie = 0; ie < num_local_elems; ++ie) {
      for (int igp = 0; igp < num_gp; ++igp) {
        for (int jgp = 0; jgp < num_gp; ++jgp) {
          int idof = ie*num_gp*num_gp + igp*num_gp + jgp;
          int gid = offset + idof;
          host_dofs(idof) = gid;
          host_lid2idx(idof, 0) = ie;
          host_lid2idx(idof, 1) = igp;
          host_lid2idx(idof, 2) = jgp;
        }
      }
    }

    dof_gids.sync_to_dev();
    lid2idx.sync_to_dev();


    add_grid(se_grid);
  }
  if (build_pt) {
    const int num_global_cols  = m_params.get<int>("number_of_global_columns");
    auto pt_grid = create_point_grid("Point Grid",num_global_cols,num_vertical_levels,m_comm);

    using namespace ShortFieldTagsNames;
    const auto units = ekat::units::Units::nondimensional();
    FieldLayout layout_mid ({LEV},{num_vertical_levels});

    auto area = pt_grid->create_geometry_data("area", pt_grid->get_2d_scalar_layout(), units);
    auto lat  = pt_grid->create_geometry_data("lat" , pt_grid->get_2d_scalar_layout(), units);
    auto lon  = pt_grid->create_geometry_data("lon" , pt_grid->get_2d_scalar_layout(), units);
    auto hyam = pt_grid->create_geometry_data("hyam", layout_mid, units);
    auto hybm = pt_grid->create_geometry_data("hybm", layout_mid, units);

    // Estimate cell area for a uniform grid by taking the surface area
    // of the earth divided by the number of columns.  Note we do this in
    // units of radians-squared.
    const Real pi        = 3.14159265358979323;
    const Real cell_area = 4.0*pi/num_global_cols;

    area.deep_copy(cell_area);
    area.sync_to_host();

    const auto nan = ekat::ScalarTraits<Real>::invalid();
    if (m_params.isParameter("latlon_filename")) {
      load_lat_lon(pt_grid);
    } else {
      lat.deep_copy(nan);
      lon.deep_copy(nan);
      lat.sync_to_host();
      lon.sync_to_host();
    }

    if (m_params.isParameter("hybrid_coefficients_filename")) {
      load_hybrid_coeffs (pt_grid);
    } else {
      hyam.deep_copy(nan);
      hybm.deep_copy(nan);
      hyam.sync_to_host();
      hybm.sync_to_host();
    }

    add_grid(pt_grid);
    this->alias_grid("Point Grid", "Physics");
  }
}

void MeshFreeGridsManager::
load_lat_lon (const nonconstgrid_ptr_type& grid) const
{
  // Get lat/lon fields
  const auto& lat = grid->get_geometry_data("lat");
  const auto& lon = grid->get_geometry_data("lon");

  using view_1d_host = AtmosphereInput::view_1d_host;
  std::map<std::string,view_1d_host> host_views = {
    { "lat", lat.template get_view<Real*,Host>() },
    { "lon", lon.template get_view<Real*,Host>() }
  };
  const auto layout = grid->get_2d_scalar_layout();
  std::map<std::string,FieldLayout> layouts = {
    { "lat", layout },
    { "lon", layout }
  };

  ekat::ParameterList reader_params;
  reader_params.set<std::vector<std::string>>("Field Names",{"lat", "lon"});
  reader_params.set("Filename",m_params.get<std::string>("latlon_filename"));

  AtmosphereInput reader(m_comm,reader_params);
  reader.init(grid,host_views,layouts);
  reader.read_variables();
  reader.finalize();

  // Make sure data is on device
  lat.sync_to_dev();
  lon.sync_to_dev();
}

void MeshFreeGridsManager::
load_hybrid_coeffs (const nonconstgrid_ptr_type& grid) const
{
  using view_1d_host = AtmosphereInput::view_1d_host;
  using namespace ShortFieldTagsNames;

  // Get hyam/hybm fields
  auto hyam = grid->get_geometry_data("hyam");
  auto hybm = grid->get_geometry_data("hybm");

  const auto nlev = grid->get_num_vertical_levels();

  std::map<std::string,view_1d_host> host_views = {
    { "hyam", hyam.template get_view<Real*,Host>() },
    { "hybm", hybm.template get_view<Real*,Host>() }
  };
  std::map<std::string,FieldLayout> layouts = {
    { "hyam", FieldLayout({LEV}, {nlev}) },
    { "hybm", FieldLayout({LEV}, {nlev}) }
  };

  ekat::ParameterList reader_params;
  reader_params.set("Filename",m_params.get<std::string>("hybrid_coefficients_filename"));
  reader_params.set<std::vector<std::string>>("Field Names",{"hyam", "hybm"});

  AtmosphereInput reader(m_comm,reader_params);
  reader.init(grid, host_views, layouts);
  reader.read_variables();
  reader.finalize();

  hyam.sync_to_dev();
  hybm.sync_to_dev();
}

std::shared_ptr<GridsManager>
create_mesh_free_grids_manager (const ekat::Comm& comm, const int num_local_elems,
                                const int num_gp, const int num_vertical_levels,
                                const int num_global_cols)
{
  ekat::ParameterList gm_params;
  gm_params.set("number_of_global_columns",num_global_cols);
  gm_params.set("number_of_local_elements",num_local_elems);
  gm_params.set("number_of_gauss_points",num_gp);
  gm_params.set("number_of_vertical_levels",num_vertical_levels);
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  return gm;
}

} // namespace scream
