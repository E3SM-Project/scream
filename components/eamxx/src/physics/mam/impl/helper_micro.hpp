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
  std::vector<view_2d> views_horiz_transpose;
  std::vector<view_2d> views_vert;
  // work arrays
  view_int_1d kupper;
  //
  view_2d pin;
  view_1d col_latitudes;
  };

  // Linoz structures to help manage all of the variables:
  struct LinozTimeState {
    // Whether the timestate has been initialized.
    // The current month
    int current_month = -1;
    // Julian Date for the beginning of the month, as defined in
    //           /src/share/util/scream_time_stamp.hpp
    // See this file for definition of Julian Date.
    Real t_beg_month;
    // Current simulation Julian Date
    Real t_now;
    // Number of days in the current month, cast as a Real
    Real days_this_month;
  }; // LinozTimeState

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


 static void perform_horizontal_interpolation( const LinozReaderParams& linoz_params)
 {
  // FIXME: get this inputs from eamxx interface.
  const auto col_latitudes = linoz_params.col_latitudes;
  const int ncol = col_latitudes.extent(0);
  const int nlev =linoz_params.views_horiz[0].extent(0);


  // We can ||ize over columns as well as over variables and bands
  const int num_vars = linoz_params.nlevs;
  LIV horiz_interp(num_vars, linoz_params.nlat, ncol);
  const auto policy_setup = ESU::get_default_team_policy(num_vars, ncol);
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


  // FIXME: Does ekat or kokkos have a call to transpose a view?
  const auto policy_transpose = ESU::get_default_team_policy(ncol*nlev,1);
  Kokkos::parallel_for("transpose_horization_data", policy_transpose,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
    const int icol = team.league_rank() / nlev;
    const int ilev = team.league_rank() % nlev;
    const auto input = linoz_params.views_horiz[0];
    const auto output = linoz_params.views_horiz_transpose[0];
    output(icol, ilev) = input(ilev, icol);
  });
  Kokkos::fence();

 }


void static
perform_vertical_interpolation_linear(
  const LinozReaderParams& linoz_params,
  const view_2d& p_tgt,
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

  // Now use the interpolation object in || over all variables.
  const int outer_iters = ncols*num_vars;
  const auto policy_interp = ESU::get_default_team_policy(outer_iters, num_vert_packs);
  Kokkos::parallel_for("spa_vert_interp_loop", policy_interp,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {

    const int icol = team.league_rank();
    const auto x1 = lev;
    const auto x2 = ekat::subview(p_tgt,icol);

    const auto y1 = ekat::subview(linoz_params.views_horiz_transpose[0],icol);
    const auto y2 = ekat::subview(linoz_params.views_vert[0],icol);
    vert_interp.lin_interp(team, x1, x2, y1, y2, icol);
  });
  Kokkos::fence();
}

// Direct port of components/eam/src/chemistry/utils/tracer_data.F90/vert_interp
static void vert_interp(int ncol,
                 int levsiz,
                 int pver,
                 const view_2d&  pin,
                 const const_view_2d&  pmid,
                 const view_2d&  datain,
                 const view_2d&  dataout,
                 //work array
                 const view_int_1d& kupper
                 ) {
    const int one = 1;
    // Initialize index array
    for (int i = 0; i < ncol; ++i) {
      kupper(i)= one;
    } // ncol

    for (int k = 0; k < pver; ++k) {
        // Top level we need to start looking is the top level for the previous k for all column points
        int kkstart = levsiz;
        for (int i = 0; i < ncol; ++i) {
            kkstart = haero::min(kkstart, kupper(i));
        }

        // Store level indices for interpolation
        for (int kk = kkstart - 1; kk < levsiz - 1; ++kk) {
            for (int i = 0; i < ncol; ++i) {
                if (pin(i, kk) < pmid(i, k) && pmid(i, k) <= pin(i, kk + 1)) {
                    kupper(i) = kk;
                }// end if
            } // end for
        } // end kk
        // Interpolate or extrapolate...
        for (int i = 0; i < ncol; ++i) {
            if (pmid(i, k) < pin(i, 0)) {
                dataout(i, k) = datain(i, 0) * pmid(i, k) / pin(i, 0);
            } else if (pmid(i, k) > pin(i, levsiz - 1)) {
                dataout(i, k) = datain(i, levsiz - 1);
            } else {
                Real dpu = pmid(i, k) - pin(i, kupper(i));
                Real dpl = pin(i, kupper(i) + 1) - pmid(i, k);
                dataout(i, k) = (datain(i, kupper[i]) * dpl + datain(i, kupper(i) + 1) * dpu) / (dpl + dpu);
            }// end if
        } // end col
    } // end k

} // vert_interp

static void perform_vertical_interpolation(const LinozReaderParams& linoz_params,
                                           const const_view_2d& p_mid)
{
  const int ncol = p_mid.extent(0);
  const int nlev = p_mid.extent(1);

  const auto kupper = linoz_params.kupper;
  const auto levs = linoz_params.io_fields[1].get_view< Real*>();
  const auto pin = linoz_params.pin;

  for (int kk = 0; kk < linoz_params.nlevs; ++kk)
  {
    const auto pin_kk = Kokkos::subview(pin,Kokkos::ALL,kk);
    Kokkos::deep_copy(pin_kk,levs(kk));
  }// i
  Kokkos::fence();

  const auto policy_interp = ESU::get_default_team_policy(1, 1);
  Kokkos::parallel_for("vertical_interpolation_linoz", policy_interp,
    KOKKOS_LAMBDA(typename LIV::MemberType const& team) {
  scream::mam_coupling::vert_interp(ncol,
              linoz_params.nlevs,
              nlev,
              pin,
              p_mid,
              linoz_params.views_horiz_transpose[0],
              linoz_params.views_vert[0],
              //work array
              kupper);
    });
    Kokkos::fence();

}//perform_vertical_interpolation

// This function is based on update_spa_timestate
void static update_linoz_timestate(const util::TimeStamp& ts,
                                   LinozTimeState& time_state,
                                   std::shared_ptr<AtmosphereInput> linoz_reader,
                                   const LinozReaderParams& linoz_params)
{
  // Now we check if we have to update the data that changes monthly
  // NOTE:  This means that SPA assumes monthly data to update.  Not
  //        any other frequency.
  const auto month = ts.get_month() - 1; // Make it 0-based
  if (month != time_state.current_month) {
    // Update the SPA time state information
    time_state.current_month = month;
    time_state.t_beg_month = util::TimeStamp({ts.get_year(),month+1,1}, {0,0,0}).frac_of_year_in_days();
    time_state.days_this_month = util::days_in_month(ts.get_year(),month+1);

    // // Copy spa_end'data into spa_beg'data, and read in the new spa_end
    // std::swap(spa_beg,spa_end);

    // Update the SPA forcing data for this month and next month
    // Start by copying next months data to this months data structure.
    // NOTE: If the timestep is bigger than monthly this could cause the wrong values
    //       to be assigned.  A timestep greater than a month is very unlikely so we
    //       will proceed.
    int next_month = (time_state.current_month + 1) % 12;
    linoz_reader->read_variables(next_month);
    perform_horizontal_interpolation(linoz_params);
  }
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
