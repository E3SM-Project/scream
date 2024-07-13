#ifndef EAMXX_MAM_HELPER_MICRO
#define EAMXX_MAM_HELPER_MICRO

#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/grid/point_grid.hpp"
#include <ekat/util/ekat_lin_interp.hpp>


namespace scream::mam_coupling {

using namespace ShortFieldTagsNames;

  struct LinozReaderParams {
  // std::shared_ptr<AtmosphereInput> linoz_reader;
  int nlevs{-1};
  int nlat{-1};
  std::vector<Field> io_fields;
  std::vector<view_2d> views_horiz;
  std::vector<view_2d> views_vert;
  };

    // using npack equal to 1.
  using LIV = ekat::LinInterp<Real,1>;
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;


  // define the different field layouts that will be used for this process
  static std::shared_ptr<AtmosphereInput>
  create_linoz_data_reader (
      const std::string& linoz_data_file,
      // std::vector<Field> io_fields,
      LinozReaderParams& linoz_params,
      const ekat::Comm& comm)
  {

     auto make_layout = [](const std::vector<int>& extents,
                        const std::vector<std::string>& names)
  {
    std::vector<FieldTag> tags(extents.size(),CMP);
    return FieldLayout(tags,extents,names);
  };

  // query the file for its resolution
  scorpio::register_file(linoz_data_file,scorpio::Read);
  const int nlevs_data = scorpio::get_dimlen(linoz_data_file,"lev");
  const int nlat_data = scorpio::get_dimlen(linoz_data_file,"lat");
  scorpio::release_file(linoz_data_file);
  std::cout << "nlevs_data: " << nlevs_data << "\n";
  std::cout << "nlat_data: " << nlat_data << "\n";
  linoz_params.nlevs=nlevs_data;
  linoz_params.nlat=nlat_data;


  // create an IO grid, with that number of cols
  // linoz files do not have number of cols,
  // I will use nlat_data instead.

  std::string name="linoz_grid";
  const auto io_grid = std::make_shared<PointGrid>(name,nlat_data,nlevs_data,comm);

  const auto nondim = ekat::units::Units::nondimensional();

  auto scalar2d_layout_linoz = make_layout({nlevs_data, nlat_data},
                                             {"lev","lat"});
  auto scalar1d_lat_layout_linoz = make_layout({nlat_data},
                                             {"lat"});
  auto scalar1d_lev_layout_linoz = make_layout({nlevs_data},
                                             {"lev"});

  Field o3_clim (FieldIdentifier("o3_clim",  scalar2d_layout_linoz,  nondim,io_grid->name()));
  Field lat (FieldIdentifier("lat",  scalar1d_lat_layout_linoz,  nondim,io_grid->name()));
  Field lev (FieldIdentifier("lev",  scalar1d_lev_layout_linoz,  nondim,io_grid->name()));
  o3_clim.allocate_view();
  lat.allocate_view();
  lev.allocate_view();
  linoz_params.io_fields.push_back(lat);
  linoz_params.io_fields.push_back(lev);
  linoz_params.io_fields.push_back(o3_clim);

  return std::make_shared<AtmosphereInput>(linoz_data_file,io_grid,linoz_params.io_fields,true);
  }


 static void perform_horizontal_interpolation( LinozReaderParams& linoz_params,
                                  const view_1d& col_latitudes)
 {
  const int ncol = col_latitudes.extent(0);

  // We can ||ize over columns as well as over variables and bands
  const int num_vars = linoz_params.nlevs;
  LIV horiz_interp(num_vars, linoz_params.nlat, ncol);
  const auto policy_setup = ESU::get_default_team_policy(num_vars, ncol);
  // Setup the linear interpolation object
  // auto& col_latitudes= col_latitudes_;
  // view_1d col_latitudes("col",ncol_);
  // Kokkos::deep_copy(col_latitudes, col_latitudes_);

  // view_1d lat("lat",nlat_data );
  // Kokkos::deep_copy(lat, lat_host);
  auto lat = linoz_params.io_fields[0].get_view<Real*>();

  Kokkos::parallel_for("spa_vert_interp_setup_loop", policy_setup,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
    // Setup
    horiz_interp.setup(team, lat, col_latitudes);
  });
  Kokkos::fence();

  // Now use the interpolation object in || over all variables.
  const int outer_iters = linoz_params.nlevs;
  const auto policy_interp = ESU::get_default_team_policy(outer_iters, ncol);
  Kokkos::parallel_for("linoz_horizongal_interp_loop", policy_interp,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
    const int kk = team.league_rank();
    const auto x1 = lat;
    const auto x2 = col_latitudes;
    auto var_org=linoz_params.io_fields[2].get_view<Real**>();
    const auto y1 = ekat::subview(var_org,kk);
    const auto y2 = ekat::subview(linoz_params.views_horiz[0],kk);
    horiz_interp.lin_interp(team, x1, x2, y1, y2, kk);
  });
  Kokkos::fence();

 }


void static
perform_vertical_interpolation(
  LinozReaderParams& linoz_params,
  view_2d& p_tgt,
  const int ncols,
  const int nlevs)
{
  const int nlevs_src = linoz_params.nlevs;
  const int nlevs_tgt = nlevs;

  LIV vert_interp(ncols,nlevs_src,nlevs_tgt);

  // We can ||ize over columns as well as over variables and bands
  const int num_vars = 1;
  const int num_vert_packs = nlevs_tgt;
  const auto policy_setup = ESU::get_default_team_policy(ncols, num_vert_packs);

  const auto lev = linoz_params.io_fields[1].get_view< Real*>();

  // Setup the linear interpolation object
  Kokkos::parallel_for("spa_vert_interp_setup_loop", policy_setup,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
    const int icol = team.league_rank();
    // Setup
    vert_interp.setup(team, lev, ekat::subview(p_tgt,icol));
  });
  Kokkos::fence();
#if 1
  // Now use the interpolation object in || over all variables.
  const int outer_iters = ncols*num_vars;
  const auto policy_interp = ESU::get_default_team_policy(outer_iters, num_vert_packs);
  Kokkos::parallel_for("spa_vert_interp_loop", policy_interp,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {

    const int icol = team.league_rank();
    const auto x1 = lev;
    const auto x2 = ekat::subview(p_tgt,icol);

    const auto y1 = ekat::subview(linoz_params.views_horiz[0],icol);
    const auto y2 = ekat::subview(linoz_params.views_vert[0],icol);
    vert_interp.lin_interp(team, x1, x2, y1, y2, icol);
  });
  Kokkos::fence();
#endif
}




#if 0
 static void vertical_interpolation( LinozReaderParams& linoz_params,
                                     const view_2d& pmid)
 {
  const int nlev = pmid.extent(1);
  const int ncol = pmid.extent(0);
  const int num_vars = linoz_params.nlevs;
  LIV horiz_interp(num_vars, nlev, ncol);
  const auto policy_setup = ESU::get_default_team_policy(num_vars, ncol);

  auto lev = linoz_params.io_fields[1].get_view<Real*>();
    // We can ||ize over columns as well as over variables and bands
  LIV ver_interp(num_vars, linoz_params.nlevs, nlev);
  const auto policy_setup = ESU::get_default_team_policy(num_vars, nlev);



 }

  static std::shared_ptr<AtmosphereInput>
  create_linoz_data_reader (
      const std::string& linoz_data_file)
  {

   using view_1d_host = typename KT::view_1d<Real>::HostMirror;
   using view_2d_host = typename KT::view_2d<Real>::HostMirror;
   using strvec_t = std::vector<std::string>;


   std::map<std::string, FieldLayout> layouts_linoz;
   // const auto& fname = m_params.get<std::string>(table_name);
   scorpio::register_file(linoz_data_file,scorpio::Read);
   const int nlevs_data = scorpio::get_dimlen(linoz_data_file,"lev");
   const int nlat_data = scorpio::get_dimlen(linoz_data_file,"lat");
   scorpio::release_file(linoz_data_file);
   std::cout << "nlevs_data" << nlevs_data << "\n";
   std::cout << "nlat_data" << nlat_data << "\n";

   auto scalar2d_layout_linoz = make_layout({nlevs_data, nlat_data},
                                             {"lev","lat"});

   auto scalar1d_lat_layout_linoz = make_layout({nlat_data},
                                             {"lat"});
   ekat::ParameterList params_Linoz;
   params_Linoz.set("Filename", linoz_data_file);
   // make a list of host views
   std::map<std::string, view_1d_host> host_views_Linoz;
   params_Linoz.set("Skip_Grid_Checks", true);
   params_Linoz.set<strvec_t>(
      "Field Names",
      {"o3_clim",
      "lat"});
   view_2d_host o3_clim_host("o3_clim_host",nlevs_data, nlat_data);
   view_1d_host lat_host("lat_host", nlat_data);
   host_views_Linoz["o3_clim"] =
      view_1d_host(o3_clim_host.data(), o3_clim_host.size());
   host_views_Linoz["lat"] =
      view_1d_host(lat_host.data(), lat_host.size());

   layouts_linoz.emplace("o3_clim", scalar2d_layout_linoz);
   layouts_linoz.emplace("lat", scalar1d_lat_layout_linoz);

   return std::make_shared<AtmosphereInput>(params_Linoz, grid_, host_views_Linoz,
                                         layouts_linoz);
  }
#endif

} // namespace scream::mam_coupling
#endif //EAMXX_MAM_HELPER_MICRO

  // scream::mam_coupling::create_linoz_data_reader(linoz_file_name,linoz_params,m_comm);

  //  auto scalar2d_layout_linoz = make_layout({nlevs_data, nlat_data},
  //                                            {"lev","lat"});

  //  auto scalar1d_lat_layout_linoz = make_layout({nlat_data},
  //                                            {"lat"});
  //  ekat::ParameterList params_Linoz;
  //  params_Linoz.set("Filename", linoz_file_name);
  //  // make a list of host views
  //  std::map<std::string, view_1d_host> host_views_Linoz;
  //  params_Linoz.set("Skip_Grid_Checks", true);
  //  params_Linoz.set<strvec_t>(
  //     "Field Names",
  //     {"o3_clim",
  //     "lat"});
  //  view_2d_host o3_clim_host("o3_clim_host",nlevs_data, nlat_data);
  //  view_1d_host lat_host("lat_host", nlat_data);
  //  host_views_Linoz["o3_clim"] =
  //     view_1d_host(o3_clim_host.data(), o3_clim_host.size());
  //  host_views_Linoz["lat"] =
  //     view_1d_host(lat_host.data(), lat_host.size());

  //  layouts_linoz.emplace("o3_clim", scalar2d_layout_linoz);
  //  layouts_linoz.emplace("lat", scalar1d_lat_layout_linoz);

  //  linoz_reader_=std::make_shared<AtmosphereInput>(params_Linoz, grid_, host_views_Linoz,
  //                                        layouts_linoz);

  //  linoz_reader_ = AtmosphereInput(params_Linoz, grid_, host_views_Linoz,
  //                                        layouts_linoz);

  //  AtmosphereInput Linoz_reader(params_Linoz, grid_, host_views_Linoz,
  //                                        layouts_linoz);
    // linoz_reader_ = scream::mam_coupling::create_linoz_data_reader(linoz_file_name);
  // linoz_reader_->finalize();
#if 0
   std::cout << "o3_clim: " << o3_clim_host(0,0) << "\n";
   std::cout << "lat: " << lat_host(0) << "\n";

   std::cout << "col_latitudes_.extents(0): " << col_latitudes_.extent(0) << "\n";
   std::cout << "ncol_: " << ncol_ << "\n";



  using LIV = ekat::LinInterp<Real,1>;
  using ExeSpace = typename KT::ExeSpace;
  using ESU = ekat::ExeSpaceUtils<ExeSpace>;

  // We can ||ize over columns as well as over variables and bands
  const int num_vars = nlevs_data;
  LIV horiz_interp(num_vars, nlat_data, ncol_);
  const auto policy_setup = ESU::get_default_team_policy(num_vars, ncol_);
  // Setup the linear interpolation object
  // auto& col_latitudes= col_latitudes_;
  col_latitudes_copy_ = view_1d("col",ncol_);
  Kokkos::deep_copy(col_latitudes_copy_, col_latitudes_);

  view_1d lat("lat",nlat_data );
  Kokkos::deep_copy(lat, lat_host);

  Kokkos::parallel_for("spa_vert_interp_setup_loop", policy_setup,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
    // Setup
    horiz_interp.setup(team, lat, col_latitudes);
  });
  Kokkos::fence();

  // Now use the interpolation object in || over all variables.
  const int outer_iters = nlevs_data;
  const auto policy_interp = ESU::get_default_team_policy(outer_iters, ncol_);
  Kokkos::parallel_for("spa_vert_interp_loop", policy_interp,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
    const int kk = team.league_rank();
    const auto x1 = lat;
    const auto x2 = col_latitudes;
    const auto y1 = ekat::subview(o3_clim_host,kk);
    const auto y2 = ekat::subview(o3_clim_test,kk);
    horiz_interp.lin_interp(team, x1, x2, y1, y2, kk);
  });
  Kokkos::fence();
  std::cout << "o3_clim_host \n";
  for (int i = 0; i < nlat_data; ++i)
  {
    std::cout <<lat(i)<<" lat " << o3_clim_host(10,i)<< ",\n";
  }

  std::cout << "o3_clim_test \n";
  for (int i = 0; i < ncol_; ++i)
  {
    std::cout <<col_latitudes(i)<<" lat " << o3_clim_test(10,i)<< ",\n";
  }
#endif
