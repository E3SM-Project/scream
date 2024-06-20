#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"
#include "physics/register_physics.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"

#include <ekat/ekat_parse_yaml_file.hpp>

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <Python.h>

#include <iomanip>

namespace scream {
TEST_CASE("ml_correction-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;
  namespace py      = pybind11;

  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  parse_yaml_file(fname, ad_params);

  const auto& ts     = ad_params.sublist("time_stepping");
  const auto  dt     = ts.get<int>("time_step");
  const auto  nsteps = ts.get<int>("number_of_steps");
  const auto  t0_str = ts.get<std::string>("run_t0");
  const auto  t0     = util::str_to_time_stamp(t0_str);
  const auto  ml     = ad_params.sublist("atmosphere_processes").sublist("MLCorrection");
  const auto  ML_model_tq_path = ml.get<std::string>("ML_model_path_tq");
  const auto  ML_model_t_only_path = ml.get<std::string>("ML_model_path_temperature");
  const auto  ML_model_uv_path = ml.get<std::string>("ML_model_path_uv");
  std::cout << "ML_model_tq_path: " << ML_model_tq_path << std::endl;
  std::cout << "ML_model_t_only_path: " << ML_model_t_only_path << std::endl;
  EKAT_ASSERT_MSG(dt > 0, "Error! Time step must be positive.\n");
  EKAT_ASSERT_MSG(!(ML_model_tq_path != "None" && ML_model_t_only_path != "None"),
                   "Error! Only one of ML_model_path_tq and ML_model_path_temperature"
                   " can be specificed. \n");  
  ekat::Comm atm_comm(MPI_COMM_WORLD);

  register_physics();
  register_mesh_free_grids_manager();

  AtmosphereDriver ad;

  ad.initialize(atm_comm, ad_params, t0);

  const auto& grid = ad.get_grids_manager()->get_grid("Physics");
  const auto& field_mgr = *ad.get_field_mgr(grid->name());

  int num_cols = grid->get_num_local_dofs();
  int num_levs = grid->get_num_vertical_levels();

  const auto &qv_field = field_mgr.get_field("qv");
  const auto &qv       = qv_field.get_view<Real **, Host>();

  for(int icol = 0; icol < num_cols; ++icol) {
    for(int jlev = 0; jlev < num_levs; ++jlev) {
      Real phase     = icol * 3.14 / 2.0 / num_cols;
      Real xval      = jlev * 3.14 / 2.0 / num_levs;
      qv(icol, jlev) = (1.0 + std::sin(xval - phase)) / 2.0;
    }
  }
  qv_field.sync_to_dev();

  Real reference = 1e-4;
  int fpe_mask = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();  // required for importing numpy
  if ( Py_IsInitialized() == 0 ) {
    py::initialize_interpreter();
  }  
  py::module sys = pybind11::module::import("sys");
  sys.attr("path").attr("insert")(1, CUSTOM_SYS_PATH);
  auto py_correction = py::module::import("test_correction");
  py::object ML_model_tq = py_correction.attr("get_ML_model")(ML_model_tq_path);
  py::object ML_model_uv = py_correction.attr("get_ML_model")(ML_model_uv_path);
  
  // CPU Test
  py::object ob1  = py_correction.attr("modify_view")(
      py::array_t<Real, py::array::c_style | py::array::forcecast>(
          num_cols * num_levs, qv.data(), py::str{}),
      num_cols, num_levs, ML_model_tq, ML_model_uv);

  // GPU Handoff test with modify
  const auto &qv_dev = qv_field.get_view<Real **, Device>();
  uintptr_t ptr = reinterpret_cast<uintptr_t>(qv_dev.data());
  std::string qv_dev_dtype = typeid(qv_dev(0, 0)).name();

  py::object test_gpu_handoff = py_correction.attr("modify_view_gpu")(
      ptr,
      qv_dev_dtype,
      num_cols, num_levs, ML_model_tq, ML_model_uv);

  // GPU handoff test xarray -- cupy -- tensorflow integration
  // TODO: where to store model for load during CI?
  py::object test_gpu_handoff_xarray = py_correction.attr("gpu_handoff_real_model")(
      ptr,
      qv_dev_dtype,
      num_cols, num_levs, ML_model_tq_path);
  
  // // Testing function for checking pointer and pybind arrays are the same
  // py::object test_ptr_usage = py_correction.attr("test_ptr")(
  //   ptr,
  //   py::array_t<Real, py::array::c_style | py::array::forcecast>(
  //         num_cols * num_levs, qv.data(), py::str{}),
  //   num_cols,
  //   num_levs); // This one should be unchanged
  
  py::gil_scoped_release no_gil;
  ekat::enable_fpes(fpe_mask);
  //Check CPU update
  REQUIRE(qv(1, 10) == reference);   // This is the one that is modified
  REQUIRE(qv(0, 10) != reference);
  
  //Check GPU update
  const auto qv_dev_h = Kokkos::create_mirror_view(qv_dev);
  Kokkos::deep_copy(qv_dev_h, qv_dev);
  REQUIRE(qv_dev_h(1, 0) == reference);   // This is the one that is modified
  REQUIRE(qv_dev_h(0, 0) != reference);
  ad.finalize();
}

}  // namespace scream
