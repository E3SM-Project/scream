#include "eamxx_ml_correction_process_interface.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"
#include "share/field/field_utils.hpp"

namespace scream {
// =========================================================================================
MLCorrection::MLCorrection(const ekat::Comm &comm,
                           const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  m_ML_model_path_tq = m_params.get<std::string>("ML_model_path_tq");
  m_ML_model_path_uv = m_params.get<std::string>("ML_model_path_uv");
  m_ML_model_path_sfc_fluxes = m_params.get<std::string>("ML_model_path_sfc_fluxes");
  m_fields_ml_output_variables = m_params.get<std::vector<std::string>>("ML_output_fields");
  m_ML_correction_unit_test = m_params.get<bool>("ML_correction_unit_test");
}

// =========================================================================================
void MLCorrection::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg / kg;
  Q.set_string("kg/kg");
  constexpr int ps = Pack::n;
  m_grid                = grids_manager->get_grid("Physics");
  const auto &grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs();  // Number of columns on this rank
  m_num_levs =
      m_grid->get_num_vertical_levels();  // Number of levels per column
  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  FieldLayout scalar2d_layout{ {COL}, {m_num_cols}};
  FieldLayout scalar3d_layout_mid{{COL, LEV}, {m_num_cols, m_num_levs}};
  FieldLayout scalar3d_layout_int{{COL, ILEV}, {m_num_cols, m_num_levs+1}};
  FieldLayout horiz_wind_layout { {COL,CMP,LEV}, {m_num_cols,2,m_num_levs} };
  if (not m_ML_correction_unit_test) {
    const auto m2 = m*m;
    const auto s2 = s*s;
    auto Wm2 = W / m / m;
    auto nondim = m/m;
    add_field<Required>("phis", scalar2d_layout, m2/s2, grid_name);
    add_field<Updated>("SW_flux_dn", scalar3d_layout_int, Wm2, grid_name, ps);
    add_field<Required>("sfc_alb_dif_vis", scalar2d_layout, nondim, grid_name);
    add_field<Updated>("sfc_flux_sw_net", scalar2d_layout, Wm2, grid_name);
    add_field<Updated>("sfc_flux_lw_dn", scalar2d_layout, Wm2, grid_name);
    m_lat  = m_grid->get_geometry_data("lat");
    m_lon  = m_grid->get_geometry_data("lon");      
  }

  /* ----------------------- WARNING --------------------------------*/
  /* The following is a HACK to get things moving, we don't want to
   * add all fields as "updated" long-term.  A separate stream of work
   * is adapting the infrastructure to allow for a generic "add_field" call
   * to be used here which we can then setup using the m_fields_ml_output_variables variable
   */
  add_field<Updated>("T_mid", scalar3d_layout_mid, K, grid_name, ps);
  add_field<Updated>("qv",    scalar3d_layout_mid, Q, grid_name, "tracers", ps);
  add_field<Updated>("horiz_winds",   horiz_wind_layout,   m/s,     grid_name, ps);
  /* ----------------------- WARNING --------------------------------*/
  add_group<Updated>("tracers", grid_name, 1, Bundling::Required);
}

void MLCorrection::initialize_impl(const RunType /* run_type */) {
  fpe_mask = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();  // required for importing numpy  
  if ( Py_IsInitialized() == 0 ) {
    pybind11::initialize_interpreter();
  }
  pybind11::module sys = pybind11::module::import("sys");
  sys.attr("path").attr("insert")(1, ML_CORRECTION_CUSTOM_PATH);
  py_correction = pybind11::module::import("ml_correction");
  py_correction.attr("initialize_logging")();
  ekat::enable_fpes(fpe_mask);
}

// =========================================================================================
void MLCorrection::run_impl(const double dt) {
  // Start timing
  auto start = std::chrono::high_resolution_clock::now();

  // use model time to infer solar zenith angle for the ML prediction
  auto current_ts = timestamp();
  std::string datetime_str = current_ts.get_date_string() + " " + current_ts.get_time_string();

  auto h_lat  = m_lat.get_view<const Real*,Host>();
  auto h_lon  = m_lon.get_view<const Real*,Host>();

  const auto& tracers = get_group_out("tracers");
  const auto& tracers_info = tracers.m_info;
  Int num_tracers = tracers_info->size();

  // Get device views instead
  const auto &qv_dev              = get_field_out("qv").get_view<Real **, Device>();
  const auto &T_mid_dev           = get_field_out("T_mid").get_view<Real **, Device>();
  const auto &phis_dev            = get_field_in("phis").get_view<const Real *, Device>();
  const auto &SW_flux_dn_dev      = get_field_out("SW_flux_dn").get_view<Real **, Device>();
  const auto &sfc_alb_dif_vis_dev = get_field_in("sfc_alb_dif_vis").get_view<const Real *, Device>();  
  const auto &sfc_flux_sw_net_dev = get_field_out("sfc_flux_sw_net").get_view<Real *, Device>();
  const auto &sfc_flux_lw_dn_dev  = get_field_out("sfc_flux_lw_dn").get_view<Real *, Device>();
  const auto &u_dev               = get_field_out("horiz_winds").get_component(0).get_view<Real **, Device>();
  const auto &v_dev               = get_field_out("horiz_winds").get_component(1).get_view<Real **, Device>();
  auto h_lat_dev  = m_lat.get_view<const Real*,Device>();
  auto h_lon_dev  = m_lon.get_view<const Real*,Device>();

  uintptr_t qv_dev_ptr = reinterpret_cast<uintptr_t>(qv_dev.data());
  std::string field_dtype = typeid(qv_dev(0, 0)).name();
  uintptr_t T_mid_dev_ptr = reinterpret_cast<uintptr_t>(T_mid_dev.data());
  uintptr_t phis_dev_ptr = reinterpret_cast<uintptr_t>(phis_dev.data());
  uintptr_t SW_flux_dn_dev_ptr = reinterpret_cast<uintptr_t>(SW_flux_dn_dev.data());
  uintptr_t sfc_alb_dif_vis_dev_ptr = reinterpret_cast<uintptr_t>(sfc_alb_dif_vis_dev.data());
  uintptr_t sfc_flux_sw_net_dev_ptr = reinterpret_cast<uintptr_t>(sfc_flux_sw_net_dev.data());
  uintptr_t sfc_flux_lw_dn_dev_ptr = reinterpret_cast<uintptr_t>(sfc_flux_lw_dn_dev.data());
  uintptr_t u_dev_ptr = reinterpret_cast<uintptr_t>(u_dev.data());
  uintptr_t v_dev_ptr = reinterpret_cast<uintptr_t>(v_dev.data());
  uintptr_t h_lat_dev_ptr = reinterpret_cast<uintptr_t>(h_lat_dev.data());
  uintptr_t h_lon_dev_ptr = reinterpret_cast<uintptr_t>(h_lon_dev.data());

  ekat::disable_all_fpes();  // required for importing numpy
  if ( Py_IsInitialized() == 0 ) {
    pybind11::initialize_interpreter();
  }
  // for qv, we need to stride across number of tracers
  pybind11::object ob1     = py_correction.attr("update_fields")(
      field_dtype, 
      qv_dev_ptr, 
      T_mid_dev_ptr,
      u_dev_ptr,
      v_dev_ptr,
      h_lat_dev_ptr,
      h_lon_dev_ptr,
      phis_dev_ptr,
      SW_flux_dn_dev_ptr,
      sfc_alb_dif_vis_dev_ptr,
      sfc_flux_sw_net_dev_ptr,
      sfc_flux_lw_dn_dev_ptr,
      m_num_cols, m_num_levs, num_tracers, dt, 
      datetime_str,
      m_ML_model_path_tq, 
      m_ML_model_path_uv, 
      m_ML_model_path_sfc_fluxes);
  pybind11::gil_scoped_release no_gil;  
  ekat::enable_fpes(fpe_mask);


  // End timing
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "MLCorrection::run_impl took " << elapsed.count() << " seconds" << std::endl;

  // Print the current call timestamp string
  std::cout << "Current timestamp: " << datetime_str << std::endl;
}

// =========================================================================================
void MLCorrection::finalize_impl() {
  // Do nothing
}
// =========================================================================================

}  // namespace scream 
