#include "dynamics/homme/homme_grids_manager.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"
#include "dynamics/homme/physics_dynamics_remapper.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"

#ifndef NDEBUG
#include "share/property_checks/field_nan_check.hpp"
#endif

#include "share/io/scorpio_input.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/remap/inverse_remapper.hpp"

// Get all Homme's compile-time dims and constants
#include "homme_dimensions.hpp"
#include "PhysicalConstants.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

namespace scream
{

HommeGridsManager::
HommeGridsManager (const ekat::Comm& comm,
                   const ekat::ParameterList& p)
 : m_comm (comm)
 , m_params(p)
{
  if (!is_parallel_inited_f90()) {
    // While we're here, we can init homme's parallel session
    auto fcomm = MPI_Comm_c2f(comm.mpi_comm());
    init_parallel_f90 (fcomm);
  }

  // This class needs Homme's context, so register as a user
  HommeContextUser::singleton().add_user();

  if (!is_params_inited_f90()) {
    // While we're here, we can init homme's parameters
    auto nlname = m_params.get<std::string>("dynamics_namelist_file_name").c_str();
    init_params_f90 (nlname);
  }

  // Check that the global number of 2d elements is no less than the number of MPI ranks
  EKAT_REQUIRE_MSG (get_homme_param<int>("nelem")>=comm.size(),
      "Error! We do not yet support running EAMxx with a number of MPI ranks\n"
      "       larger than the number of 2d elements in Homme.\n"
      "  - num MPI ranks: " + std::to_string(comm.size()) + "\n"
      "  - num 2d elems : " + std::to_string(get_homme_param<int>("nelem")) + "\n");

  // Create the grid integer codes map (i.e., int->string
  build_pg_codes ();
}

HommeGridsManager::
~HommeGridsManager () {
  // Cleanup the grids stuff
  finalize_geometry_f90 ();

  // This class is done with Homme. Remove from its users list
  HommeContextUser::singleton().remove_user();
}

HommeGridsManager::remapper_ptr_type
HommeGridsManager::do_create_remapper (const grid_ptr_type from_grid,
                                       const grid_ptr_type to_grid) const {
  const auto from = from_grid->name();
  const auto to   = to_grid->name();

  EKAT_REQUIRE_MSG ( from=="Dynamics" || to=="Dynamics",
    "Error! Either source or target grid must be 'Dynamics'.\n");

  const bool p2d = to=="Dynamics";

  if (from=="Physics GLL" || to=="Physics GLL") {
    using PDR = PhysicsDynamicsRemapper;

    auto dyn_grid = get_grid("Dynamics");
    auto phys_grid = get_grid("Physics GLL");

    auto pd_remapper = std::make_shared<PDR>(phys_grid,dyn_grid);
    if (p2d) {
      return pd_remapper;
    } else {
      return std::make_shared<InverseRemapper>(pd_remapper);
    }
  } else {
    ekat::error::runtime_abort("Error! P-D remapping only implemented for 'Physics GLL' phys grid.\n");
  }
  return nullptr;
}

void HommeGridsManager::
build_grids ()
{
  // Get the physics grid specs
  const ci_string pg_type      = m_params.get<std::string>("physics_grid_type");
  const ci_string pg_rebalance = m_params.get<std::string>("physics_grid_rebalance","None");

  // Get the physics grid code
  std::vector<int> pg_codes {
    m_pg_codes["GLL"]["None"],  // We always need this to read/write dyn grid stuff
    m_pg_codes[pg_type][pg_rebalance]
  };
  // In case the two pg codes are the same...
  auto it = std::unique(pg_codes.begin(),pg_codes.end());
  const int* codes_ptr = pg_codes.data();
  init_grids_f90 (codes_ptr,std::distance(pg_codes.begin(),it));
  
  // We know we need the dyn grid, so build it
  build_dynamics_grid ();

  // Also the GLL grid with no rebalance is needed for sure
  build_physics_grid("GLL","None");

  // If (pg type,rebalance) is (GLL,None), this will be a no op
  build_physics_grid(pg_type,pg_rebalance);

  // Make "Physics" be an alias to whatever the pair (pg_type,rebalance) refers to
  std::string pg_name = "Physics " + pg_type;
  if (pg_rebalance!="None") {
    pg_name += " " + pg_rebalance;
  }

  this->alias_grid(pg_name,"Physics");

  // Clean up temporaries used during grid initialization
  cleanup_grid_init_data_f90 ();
}

void HommeGridsManager::build_dynamics_grid () {
  const std::string name = "Dynamics";
  if (has_grid(name)) {
    return;
  }

  using gid_t = AbstractGrid::gid_type;

  // Get dimensions and create "empty" grid
  const int nlelem = get_num_local_elems_f90();
  const int nlev   = get_nlev_f90();

  auto dyn_grid = std::make_shared<SEGrid>("Dynamics",nlelem,HOMMEXX_NP,nlev,m_comm);
  dyn_grid->setSelfPointer(dyn_grid);

  const auto layout2d = dyn_grid->get_2d_scalar_layout();
  const auto rad = ekat::units::Units::nondimensional();

  // Filling the cg/dg gids, elgpgp, coords, lat/lon views
  auto dg_dofs = dyn_grid->get_dofs_gids();
  auto cg_dofs = dyn_grid->get_cg_dofs_gids();
  auto elgpgp  = dyn_grid->get_lid_to_idx_map();
  auto lat     = dyn_grid->create_geometry_data("lat",layout2d,rad);
  auto lon     = dyn_grid->create_geometry_data("lon",layout2d,rad);

  auto dg_dofs_h = dg_dofs.get_view<gid_t*,Host>();
  auto cg_dofs_h = cg_dofs.get_view<gid_t*,Host>();
  auto elgpgp_h  = elgpgp.get_view<int**,Host>();
  auto lat_h     = lat.get_view<Real***,Host>();
  auto lon_h     = lon.get_view<Real***,Host>();

  // Get (ie,igp,jgp,gid) data for each dof
  get_dyn_grid_data_f90 (dg_dofs_h.data(),cg_dofs_h.data(),elgpgp_h.data(), lat_h.data(), lon_h.data());

  dg_dofs.sync_to_dev();
  cg_dofs.sync_to_dev();
  elgpgp.sync_to_dev();
  lat.sync_to_dev();
  lon.sync_to_dev();

#ifndef NDEBUG
  auto lat_check = std::make_shared<FieldNaNCheck>(lat,dyn_grid)->check();
  EKAT_REQUIRE_MSG (lat_check.result==CheckResult::Pass,
      "ERROR! NaN values detected in latitude field.\n" + lat_check.msg);
  auto lon_check = std::make_shared<FieldNaNCheck>(lon,dyn_grid)->check();
  EKAT_REQUIRE_MSG (lon_check.result==CheckResult::Pass,
      "ERROR! NaN values detected in longitude field.\n" + lon_check.msg);
#endif

  init_homme_vcoord(dyn_grid);

  add_grid(dyn_grid);
}

void HommeGridsManager::
build_physics_grid (const ci_string& type, const ci_string& rebalance) {
  std::string name = "Physics " + type;
  if (rebalance != "None") {
    name += " " + rebalance;
  }

  // Build only if not built yet
  if (has_grid(name)) {
    return;
  }

  // Get the grid pg_type
  const int pg_code = m_pg_codes.at(type).at(rebalance);

  // Get dimensions and create "empty" grid
  const int nlev  = get_nlev_f90();
  const int nlcols = get_num_local_columns_f90 (pg_code % 10);

  auto phys_grid = std::make_shared<PointGrid>(name,nlcols,nlev,m_comm);
  phys_grid->setSelfPointer(phys_grid);

  // Create the gids, coords, area views
  using namespace ShortFieldTagsNames;
  const auto layout2d = phys_grid->get_2d_scalar_layout();
  const auto nondim = ekat::units::Units::nondimensional();

  auto dofs = phys_grid->get_dofs_gids();
  auto lat  = phys_grid->create_geometry_data("lat",layout2d,nondim);
  auto lon  = phys_grid->create_geometry_data("lon",layout2d,nondim);
  auto area = phys_grid->create_geometry_data("area",layout2d,nondim);

  using gid_t = AbstractGrid::gid_type;

  auto dofs_h = dofs.get_view<gid_t*,Host>();
  auto lat_h  = lat.get_view<Real*,Host>();
  auto lon_h  = lon.get_view<Real*,Host>();
  auto area_h = lon.get_view<Real*,Host>();

  // Get all specs of phys grid cols (gids, coords, area)
  get_phys_grid_data_f90 (pg_code, dofs_h.data(), lat_h.data(), lon_h.data(), area_h.data());

  dofs.sync_to_dev();
  lat.sync_to_dev();
  lon.sync_to_dev();
  area.sync_to_dev();

#ifndef NDEBUG
  auto lat_check = std::make_shared<FieldNaNCheck>(lat,phys_grid)->check();
  EKAT_REQUIRE_MSG (lat_check.result==CheckResult::Pass,
      "ERROR! NaN values detected in latitude field.\n" + lat_check.msg);
  auto lon_check = std::make_shared<FieldNaNCheck>(lon,phys_grid)->check();
  EKAT_REQUIRE_MSG (lon_check.result==CheckResult::Pass,
      "ERROR! NaN values detected in longitude field.\n" + lon_check.msg);
  auto area_check = std::make_shared<FieldNaNCheck>(area,phys_grid)->check();
  EKAT_REQUIRE_MSG (area_check.result==CheckResult::Pass,
      "ERROR! NaN values detected in area field.\n" + area_check.msg);
#endif

  // Copy hybrid coords from dyn grid into phys grid
  auto hyai = get_grid("Dynamics")->get_geometry_data("hyai");
  auto hybi = get_grid("Dynamics")->get_geometry_data("hybi");
  auto hyam = get_grid("Dynamics")->get_geometry_data("hyam");
  auto hybm = get_grid("Dynamics")->get_geometry_data("hybm");
  phys_grid->set_geometry_data(hyai);
  phys_grid->set_geometry_data(hybi);
  phys_grid->set_geometry_data(hyam);
  phys_grid->set_geometry_data(hybm);

  add_grid(phys_grid);
}

void HommeGridsManager::init_homme_vcoord (const nonconstgrid_ptr_type& dyn_grid) {
  using view_1d_host = AtmosphereInput::view_1d_host; 
  using vos_t = std::vector<std::string>;
  using namespace ShortFieldTagsNames;

  const int nlev_int = HOMMEXX_NUM_INTERFACE_LEV;
  const int nlev_mid = HOMMEXX_NUM_PHYSICAL_LEV;

  // Read vcoords into host views
  ekat::ParameterList vcoord_reader_pl;
  vcoord_reader_pl.set("Filename",m_params.get<std::string>("vertical_coordinate_filename"));
  vcoord_reader_pl.set<vos_t>("Field Names",{"hyai","hybi","hyam","hybm"});

  FieldLayout layout_mid ({ LEV},{nlev_mid});
  FieldLayout layout_int ({ILEV},{nlev_int});
  constexpr auto nondim = ekat::units::Units::nondimensional();

  auto hyai = dyn_grid->create_geometry_data("hyai",layout_int,nondim);
  auto hybi = dyn_grid->create_geometry_data("hybi",layout_int,nondim);
  auto hyam = dyn_grid->create_geometry_data("hyam",layout_mid,nondim);
  auto hybm = dyn_grid->create_geometry_data("hybm",layout_mid,nondim);

  std::map<std::string,view_1d_host> host_views = {
    { "hyai", hyai.get_view<Real*,Host>() },
    { "hybi", hybi.get_view<Real*,Host>() },
    { "hyam", hyam.get_view<Real*,Host>() },
    { "hybm", hybm.get_view<Real*,Host>() }
  };
  std::map<std::string,FieldLayout> layouts = {
    { "hyai", layout_int },
    { "hybi", layout_int },
    { "hyam", layout_mid },
    { "hybm", layout_mid }
  };

  AtmosphereInput vcoord_reader(m_comm,vcoord_reader_pl);
  vcoord_reader.init(dyn_grid,host_views,layouts);
  vcoord_reader.read_variables();
  vcoord_reader.finalize();

  // Sync to device
  hyai.sync_to_dev();
  hybi.sync_to_dev();
  hyam.sync_to_dev();
  hybm.sync_to_dev();
  
  // Pass host views data to hvcoord init function
  const auto ps0 = Homme::PhysicalConstants::p0;

  // Set vcoords in f90
  prim_set_hvcoords_f90 (ps0,
                         host_views["hyai"].data(),
                         host_views["hybi"].data(),
                         host_views["hyam"].data(),
                         host_views["hybm"].data());
}

void HommeGridsManager::
build_pg_codes () {

  // Codes for the physics grids supported
  // Rules: a code has two digits, X and Y
  //  - X denotes the column rebalance choice:
  //    - 0: no rebalance
  //    - 1: twin columns rebalance
  //  - Y denotes the type of physics grid:
  //    - 0: GLL grid
  //    - N: FV phys grid, with NxN points

  m_pg_codes["GLL"]["None"] =  0;
  m_pg_codes["GLL"]["Twin"] = 10;
  m_pg_codes["PG2"]["None"] =  2;
  m_pg_codes["PG2"]["Twin"] = 12;
} 

} // namespace scream
