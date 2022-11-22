#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/remap/do_nothing_remapper.hpp"
#include "share/io/scorpio_input.hpp"

#include "physics/share/physics_constants.hpp"

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

    // Set up the degrees of freedom.
    SEGrid::dofs_list_type dofs("", num_local_elems*num_gp*num_gp);
    SEGrid::lid_to_idx_map_type dofs_map("", num_local_elems*num_gp*num_gp, 3);

    auto host_dofs = Kokkos::create_mirror_view(dofs);
    auto host_dofs_map = Kokkos::create_mirror_view(dofs_map);

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
          host_dofs_map(idof, 0) = ie;
          host_dofs_map(idof, 1) = igp;
          host_dofs_map(idof, 2) = jgp;
        }
      }
    }

    // Move the data to the device and set the DOFs.
    Kokkos::deep_copy(dofs, host_dofs);
    Kokkos::deep_copy(dofs_map, host_dofs_map);

    // Create the grid, and set the dofs
    std::shared_ptr<SEGrid> se_grid;
    se_grid = std::make_shared<SEGrid>("SE Grid",num_local_elems,num_gp,num_vertical_levels,m_comm);
    se_grid->setSelfPointer(se_grid);

    se_grid->set_dofs(dofs);
    se_grid->set_lid_to_idx_map(dofs_map);

    add_grid(se_grid);
  }
  if (build_pt) {
    const int num_global_cols  = m_params.get<int>("number_of_global_columns");
    auto pt_grid = create_point_grid("Point Grid",num_global_cols,num_vertical_levels,m_comm);

    using geo_view_type = AbstractGrid::geo_view_type;
    using PC            = scream::physics::Constants<Real>;

    const int num_local_cols = pt_grid->get_num_local_dofs();

    // Estimate cell area for a uniform grid by taking the surface area
    // of the earth divided by the number of columns.  Note we do this in
    // units of radians-squared.
    geo_view_type area("area", num_local_cols);
    const Real pi        = PC::Pi;
    const Real cell_area = 4.0*pi/num_local_cols;
    Kokkos::deep_copy(area, cell_area);
    pt_grid->set_geometry_data("area", area);

    // Load lat/lon if latlon_filename param is given.
    if (m_params.isParameter("latlon_filename")) {
      load_lat_lon(pt_grid);
    }

    // Load hyam/hybm if hybrid_coefficients_filename param is given.
    if (m_params.isParameter("vertical_coordinate_filename")) {
      load_vertical_coordinates(pt_grid);
    }

    // Load topography from file if topography_filename
    // param is given.
    // Create topography geometry data (filled with NaNs)
    // if create_topography_without_file=true (ex. used
    // for diagnostic unit tests where a random topo is needed).
    if (m_params.isParameter("topography_filename")) {
      EKAT_REQUIRE_MSG(not m_params.get<bool>("create_topography_without_file", false),
                      "Error! Setting parameter topography_filename and "
                      "create_topography_without_file=true is not allowed.\n");
      load_topography(pt_grid);
    } else if (m_params.get<bool>("create_topography_without_file", false)) {
      geo_view_type topo("topo",num_local_cols);
      Kokkos::deep_copy(topo,ekat::ScalarTraits<Real>::invalid());
      pt_grid->set_geometry_data("topo",topo);
    }

    add_grid(pt_grid);
    this->alias_grid("Point Grid", "Physics");
  }
}

namespace {

// Helper function for creating mirror device view used to set geometry data.
const AbstractGrid::geo_view_type cmvc(const AbstractGrid::geo_view_h_type host_view) {
  return Kokkos::create_mirror_view_and_copy(DefaultDevice(), host_view);
}

}

void MeshFreeGridsManager::
load_lat_lon (const nonconstgrid_ptr_type& grid) const
{
  const int num_local_cols = grid->get_num_local_dofs();

  using geo_view_host = AbstractGrid::geo_view_type::HostMirror;

  // Create host mirrors for reading in data
  std::map<std::string,geo_view_host> host_views = {
    { "lat", geo_view_host("lat",num_local_cols) },
    { "lon", geo_view_host("lon",num_local_cols) }
  };

  // Store view layouts
  std::map<std::string,FieldLayout> layouts = {
    { "lat", grid->get_2d_scalar_layout() },
    { "lon", grid->get_2d_scalar_layout() }
  };

  // Read lat/lon into host views
  ekat::ParameterList lat_lon_reader_pl;
  lat_lon_reader_pl.set("Filename",m_params.get<std::string>("latlon_filename"));
  lat_lon_reader_pl.set<std::vector<std::string>>("Field Names",{"lat","lon"});

  AtmosphereInput lat_lon_reader(m_comm, lat_lon_reader_pl);
  lat_lon_reader.init(grid, host_views, layouts);
  lat_lon_reader.read_variables();
  lat_lon_reader.finalize();

  // Set geometry data
  grid->set_geometry_data("lat", cmvc(host_views["lat"]));
  grid->set_geometry_data("lon", cmvc(host_views["lon"]));
}

void MeshFreeGridsManager::
load_vertical_coordinates (const nonconstgrid_ptr_type& grid) const
{
  const int num_vertical_levels = grid->get_num_vertical_levels();

  using geo_view_host = AbstractGrid::geo_view_type::HostMirror;

  // Create host mirrors for reading in data
  std::map<std::string,geo_view_host> host_views = {
    { "hyam", geo_view_host("hyam",num_vertical_levels) },
    { "hybm", geo_view_host("hybm",num_vertical_levels) }
  };

  // Store view layouts
  using namespace ShortFieldTagsNames;
  std::map<std::string,FieldLayout> layouts = {
    { "hyam", FieldLayout({LEV},{num_vertical_levels}) },
    { "hybm", FieldLayout({LEV},{num_vertical_levels}) }
  };

  // Read hyam/hybm into host views
  ekat::ParameterList vcoord_reader_pl;
  vcoord_reader_pl.set("Filename",m_params.get<std::string>("vertical_coordinate_filename"));
  vcoord_reader_pl.set<std::vector<std::string>>("Field Names",{"hyam","hybm"});

  AtmosphereInput vcoord_reader(m_comm,vcoord_reader_pl);
  vcoord_reader.init(grid, host_views, layouts);
  vcoord_reader.read_variables();
  vcoord_reader.finalize();

  // Set to geometry data
  grid->set_geometry_data("hyam", cmvc(host_views["hyam"]));
  grid->set_geometry_data("hybm", cmvc(host_views["hybm"]));
}

void MeshFreeGridsManager::
load_topography (const nonconstgrid_ptr_type& grid) const
{
  const int num_local_cols = grid->get_num_local_dofs();

  using geo_view_host = AbstractGrid::geo_view_type::HostMirror;

  // TODO: Use topo files instead of IC files. This hack is currently needed since
  //       the topo files use uppercase PHIS, whereas IC files use lower case.
  //       Right now the issue is that not all CIME compsets have equiv. topo
  //       files. Ex. aquaplanet file is needed with phis=0.
#if 0
    std::string topo_name = "PHIS";
#else
    std::string topo_name = "phis";
#endif

  // Create host mirrors for reading in data
  std::map<std::string,geo_view_host> host_views = {
    { topo_name, geo_view_host("topo",num_local_cols) }
  };

  // Store view layouts
  std::map<std::string,FieldLayout> layouts = {
    { topo_name, grid->get_2d_scalar_layout() }
  };

  // Read topography into host views
  ekat::ParameterList topo_reader_pl;
  topo_reader_pl.set("Filename",m_params.get<std::string>("topography_filename"));
  topo_reader_pl.set<std::vector<std::string>>("Field Names",{topo_name});

  AtmosphereInput topo_reader(m_comm, topo_reader_pl);
  topo_reader.init(grid, host_views, layouts);
  topo_reader.read_variables();
  topo_reader.finalize();

  // Set geometry data
  grid->set_geometry_data("topo",cmvc(host_views[topo_name]));
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
