#include "catch2/catch.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "physics/spa/spa_functions.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace {

using namespace scream;
using namespace spa;

template <typename S>
using view_1d = typename KokkosTypes<DefaultDevice>::template view_1d<S>;
template <typename S>
using view_1d_host = typename view_1d<S>::HostMirror; //KokkosTypes<HostDevice>::template view_1d<S>;

using SPAFunc  = spa::SPAFunctions<Real, DefaultDevice>;
using C = scream::physics::Constants<Real>;
constexpr auto test_tol = C::macheps*10;
//=====================================================================//
//  Helper function to check if two interp structures match
//=====================================================================//
void interp_is_equal(const SPAFunc::SPAHorizInterp& spa_horiz_interp_A, const SPAFunc::SPAHorizInterp& spa_horiz_interp_B)
{
  REQUIRE(spa_horiz_interp_A.length == spa_horiz_interp_B.length);
  REQUIRE(spa_horiz_interp_A.num_unique_cols == spa_horiz_interp_B.num_unique_cols);
  // Two remappers may have the same data, but organized differently. So we simply check that any triplet
  //   [src_col, tgt_col, weight] can be found in both A and B
  // To make sure we don't match two from A to one from B we prune any matches.
  std::vector<int> B_indices(spa_horiz_interp_B.length);
  std::iota(B_indices.begin(),B_indices.end(),0);
  for (int ii=0;ii<spa_horiz_interp_A.length;ii++) {
    bool found = false;
    for (int jj=0;jj<B_indices.size();jj++) {
      int idx = B_indices[jj];
      if (
           spa_horiz_interp_A.target_grid_loc(ii) == spa_horiz_interp_B.target_grid_loc(idx) &&
           spa_horiz_interp_A.source_grid_loc(ii) == spa_horiz_interp_B.source_grid_loc(idx) &&
           std::abs(spa_horiz_interp_A.weights(ii) - spa_horiz_interp_B.weights(idx)) < test_tol
         ) {
        found = true;
        B_indices.erase(B_indices.begin()+jj);
        break;
      }
    }
    REQUIRE(found);
  }
}

TEST_CASE("spa_read_remap_data","spa")
{
  //======================= Set up the mpi communicator and init the pio subsystem =======================// 
  ekat::Comm spa_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  MPI_Fint fcomm = MPI_Comm_c2f(spa_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  //======================= Establish the SPA function object =======================//
  using gid_type = SPAFunc::gid_type;
  SPAFunc::SPAHorizInterp spa_horiz_interp;
  spa_horiz_interp.m_comm = spa_comm;

  //======================= Hardcoded location for spa data file for testing. =======================//
  const std::string remap_file_name = SCREAM_DATA_DIR "/init/spa_data_for_testing_20220801.nc";

  //======================= Setup test size =======================//
  Int tgt_grid_ncols_total = 48;
  Int src_grid_ncols = 20;

  //======================= Break the test set of columns into local degrees of freedom per mpi rank =======================//
  auto comm_size = spa_comm.size();
  auto comm_rank = spa_comm.rank();
  int tgt_grid_ncols = tgt_grid_ncols_total/comm_size + (comm_rank < tgt_grid_ncols_total%comm_size ? 1 : 0);
  view_1d<gid_type> dofs_gids("",tgt_grid_ncols);
  gid_type min_dof = 0;  // We will set up the dof's to start from 0
  Kokkos::parallel_for("", tgt_grid_ncols, KOKKOS_LAMBDA(const int& ii) {
    dofs_gids(ii) = min_dof + static_cast<gid_type>(comm_rank + ii*comm_size);
  });
  // Build a map of the dof and local index on this rank, needed for testing
  std::map<int,int> local_dof_map;
  for (int ii=0;ii<tgt_grid_ncols;ii++) {
    local_dof_map.emplace(static_cast<int>(dofs_gids(ii)),ii);
  }

  //======================= Make sure that the total set of columns has been completely broken up. =======================//
  Int test_total_ncols = 0;
  spa_comm.all_reduce(&tgt_grid_ncols,&test_total_ncols,1,MPI_SUM);
  REQUIRE(test_total_ncols == tgt_grid_ncols_total);

  //======================= Test get_remap_indices =======================//
  std::vector<int> seg_dof, seg_start, seg_length;
  SPAFunc::get_remap_indices(spa_comm,remap_file_name,dofs_gids,seg_dof,seg_start,seg_length);
  // Check that all segments on this rank match the DOF's on this rank:
  for (int ii=0;ii<seg_dof.size();ii++) {
    bool match = false;
    int tgt_col = seg_dof[ii];
    for (int jj=0;jj<tgt_grid_ncols;jj++) {
      if (dofs_gids(jj)==tgt_col) {
        match = true; 
        break;
      }
    }
    REQUIRE(match);
  }

  //======================= Test get_remap_weights_from_file =======================//
  SPAFunc::get_remap_weights_from_file(remap_file_name,spa_horiz_interp,seg_dof,seg_start,seg_length);
  {
    // Check that the total size of data matches the amount that is in the file
    int total_length;
    spa_comm.all_reduce(&spa_horiz_interp.length,&total_length,1,MPI_SUM);
    scorpio::register_file(remap_file_name,scorpio::Read);
    int true_length = scorpio::get_dimlen_c2f(remap_file_name.c_str(),"n_s");
    scorpio::eam_pio_closefile(remap_file_name);
    REQUIRE(total_length==true_length);
  }
  // Check that the total weights all add up to one, for each DOF
  {
    // Create a view of weights for all the dof's on this rank
    view_1d<Real> weights("",tgt_grid_ncols);
    Kokkos::deep_copy(weights,0.0);
    for (int ii=0;ii<spa_horiz_interp.length;ii++) {
      const int dof_i = spa_horiz_interp.target_grid_loc(ii);
      const int idx   = local_dof_map[dof_i];
      weights(idx) += spa_horiz_interp.weights(ii);
    }
    for (int ii=0;ii<tgt_grid_ncols;ii++) {
      EKAT_REQUIRE_MSG(std::abs(weights(ii)-1.0)<test_tol,"Error: spa remap weights don't sum to 1 for dof = " 
                             + std::to_string(dofs_gids(ii)) + ", weight = " + std::to_string(weights(ii)) +".");
    }
  }
  // Check that the source -> target setup matches the test.  To do this we create the SpaHorizInterp we would
  // expect and check that it matches the one from get_remap_weights_from_file.
  // TODO: It would be great to create function that just constructs the remap file from scratch, that way
  //       we can be confident of the structure.  For now, the file is static following the convention
  //       below.
  {
    int length = 0;
    std::vector<int> src_grid_loc, tgt_grid_loc;
    for (int ii=0;ii<tgt_grid_ncols;ii++) {
      const int tgt_col = dofs_gids(ii);
      const int src_max = std::min(tgt_col+1,src_grid_ncols);
      for (int jj=0;jj<src_max;jj++) {
        const int idx = length;
        src_grid_loc.push_back(jj);
        tgt_grid_loc.push_back(tgt_col);
        length += 1;
      }
    }
    SPAFunc::SPAHorizInterp spa_horiz_interp_test(length);
    Kokkos::parallel_for("", length, [&] (const int& ii) {
      spa_horiz_interp_test.source_grid_loc(ii) = src_grid_loc[ii];
      spa_horiz_interp_test.target_grid_loc(ii) = tgt_grid_loc[ii];
      const int src_max = std::min(tgt_grid_loc[ii]+1, src_grid_ncols);
      spa_horiz_interp_test.weights(ii)         = (1.0+src_grid_loc[ii]) / ( (src_max) * (src_max+1.0) / 2.0 );
    });
    spa_horiz_interp_test.set_unique_cols();

    interp_is_equal(spa_horiz_interp, spa_horiz_interp_test);    
  }
  
  // All Done 
  scorpio::eam_pio_finalize();
}

} //namespace
