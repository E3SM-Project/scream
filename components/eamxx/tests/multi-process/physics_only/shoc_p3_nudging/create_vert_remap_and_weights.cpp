#include <catch2/catch.hpp>
#include "share/io/scream_output_manager.hpp"
#include "share/io/scream_scorpio_interface.hpp"

namespace {

using namespace scream;

void create_vert_remap() {
  // Simple function to create a 1D remap column to test nudging w/ remapped data

  int nlevs = 5*SCREAM_PACK_SIZE+1;
  std::vector<std::int64_t> dofs_levs(nlevs);
  std::iota(dofs_levs.begin(),dofs_levs.end(),0);
  std::vector<Real> p_tgt;
  Real p_top=0, p_bot=102500;
  Real dp = (p_bot - p_top) / (nlevs-1);
  for (int ii=0; ii<nlevs; ++ii) {
    Real p_loc = p_top + dp*ii;
    p_tgt.push_back(p_loc);
  }

  std::string remap_filename = "vertical_remap.nc";

  scorpio::register_file(remap_filename, scorpio::FileMode::Write);
  scorpio::define_dim(remap_filename,"lev",nlevs);
  scorpio::define_var(remap_filename,"p_levs",{"lev"},"real");
  scorpio::enddef(remap_filename);
  scorpio::write_var(remap_filename,"p_levs",p_tgt.data());
  scorpio::release_file(remap_filename);
}

void create_nudging_weights_ncfile(int ntimes, int ncols, int nlevs, const std::string& filename)
{
  // Simple function to create a 1D remap column to test nudging w/ remapped data
  Real plev[nlevs];
  Real p_top=0, p_bot=102500;
  Real dp = (p_bot - p_top) / (nlevs-1);
  for (int ilev=0; ilev<nlevs; ++ilev) {
    plev[ilev] = p_top + dp*ilev;
  }  

  scorpio::register_file(filename, scorpio::FileMode::Write);
  scorpio::define_dim(filename,"ncol", ncols);
  scorpio::define_dim(filename,"lev",  nlevs);
  scorpio::define_time(filename,"nsteps");
  scorpio::define_var(filename,"nudging_weights",{"ncol", "lev"},"real",true);
  scorpio::enddef(filename);

  Real weights[ncols][nlevs];
  for (auto itime=0; itime<ntimes; ++itime) {
    for (auto ilev=0; ilev<nlevs; ++ilev) {
      for (auto icol=0; icol<ncols; ++icol) {
        if (plev[ilev] <= 1.0e5 && plev[ilev] >= 8.0e4) {
           weights[icol][ilev] = 1.;
        } else {
           weights[icol][ilev] = 0.;
        }
      }
    }
    scorpio::update_time(filename,itime);
    scorpio::write_var(filename,"nudging_weights",&weights[0][0]);
  }

  scorpio::release_file(filename);
}

TEST_CASE("create_vert_remap_and_weights","create_vert_remap_and_weights")
{
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_pio_subsystem(comm);

  create_vert_remap();
  create_nudging_weights_ncfile(1, 218, 72, "nudging_weights.nc");
  create_nudging_weights_ncfile(1, 218, 128, "nudging_weights_L128.nc");

  scorpio::finalize_pio_subsystem();
}
} // end namespace
