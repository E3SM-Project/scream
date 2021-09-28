#ifndef SPA_FUNCTIONS_HPP
#define SPA_FUNCTIONS_HPP

#include "share/scream_types.hpp"
#include "share/util/scream_time_stamp.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"
#include "ekat/mpi/ekat_comm.hpp"

namespace scream {
namespace spa {

template <typename ScalarT, typename DeviceT>
struct SPAFunctions
{

  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S>
  using BigPack = ekat::Pack<S,SCREAM_PACK_SIZE>;
  template <typename S>
  using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  using WorkspaceManager = typename ekat::WorkspaceManager<Spack, Device>;
  using Workspace        = typename WorkspaceManager::Workspace;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using view_3d = typename KT::template view_3d<S>;
  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  /* ------------------------------------------------------------------------------------------- */
  // SPA structures to help manage all of the variables:
  struct SPATimeState {
    SPATimeState() = default;
    // The current month
    Int current_month;
    // Julian Date for the beginning of the month, as defined in
    //           /src/share/util/scream_time_stamp.hpp
    // See this file for definition of Julian Date.
    Real t_beg_month;
    // Current simulation Julian Date
    Real t_now;
    // Number of days in the current month, cast as a Real
    Real days_this_month;
  }; // SPATimeState

  struct SPAPressureState {
    SPAPressureState() = default;
    // Number of horizontal columns and vertical levels for the data
    Int ncols;
    Int nlevs;
    // Surface pressure for data at the beginning of the month
    view_1d<Real> ps_this_month;
    // Surface pressure for data at the beginning of next month
    view_1d<Real> ps_next_month;
    // Hybrid coordinate values
    view_1d<const Spack> hyam, hybm;
    // Current simulation pressure levels
    view_2d<const Spack> pmid;
  }; // SPAPressureState

  struct SPAData {
    SPAData() = default;
    // Number of horizontal columns and vertical levels for the data
    Int ncols;
    Int nlevs;
    // CCN3: CCN concentration at S=0.1%, units = #/cm3, dimensions = (ncol,nlev)
    view_2d<Spack> CCN3;
    // AER_G_SW: unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_G_SW;
    // AER_SSA_SW: unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_SSA_SW;
    // AER_TAU_SW: unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_TAU_SW;
    // AER_TAU_LW: unit = #/cm3, dimensions = (ncol,nswband=16,nlev)
    view_3d<Spack> AER_TAU_LW;
  }; // SPAPrescribedAerosolData

  struct SPAOutput {
    SPAOutput() = default;
    // Number of horizontal columns and vertical levels for the data
    Int ncols;
    Int nlevs;
    // CCN3: CCN concentration at S=0.1%, units = #/cm3, dimensions = (ncol,nlev)
    view_2d<Spack> CCN3;
    // AER_G_SW - 14 bands
    // AER_G_SW: unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_G_SW;
    // AER_SSA_SW: unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_SSA_SW;
    // AER_TAU_SW: unit = #/cm3, dimensions = (ncol,nswband=14,nlev)
    view_3d<Spack> AER_TAU_SW;
    // AER_TAU_LW: unit = #/cm3, dimensions = (ncol,nswband=16,nlev)
    view_3d<Spack> AER_TAU_LW;
  }; // SPAPrescribedAerosolData

  struct SPAInterp {
    SPAInterp() = default;
    // Length of remap data
    Int length;
    // 1D index of weights. Needs decoder indexing.
    view_1d<Real> weights;
    // Index of Source Grid Column.
    view_1d<Int> src_grid_loc;
    // Index of Destination Grid Column.
    view_1d<Int> dst_grid_loc;
    // Number of columns and levels on source grid
    Int src_grid_ncols, src_grid_nlevs;
  }; // SPAInterp
  /* ------------------------------------------------------------------------------------------- */
  // SPA routines
  static void spa_main(
    const SPATimeState& time_state,
    const SPAPressureState& pressure_state,
    const SPAData&   data_beg,
    const SPAData&   data_end,
    const SPAOutput& data_out,
    const Int ncols_scream,
    const Int nlevs_scream,
    const Int nswbands,
    const Int nlwbands);

  static void read_remap_weights_from_file(
    const ekat::Comm& comm,
    const std::string& remap_file,
    const SPAInterp& horiz_weights 
  );

  static void spa_update_monthly_data(
    const ekat::Comm& comm,
    const std::string& spa_data_file,
    const util::TimeStamp& ts,
    const SPAInterp& horiz_weights, 
          SPATimeState& time_state,
          SPAPressureState& pressure_state,
          SPAData&   data_beg,
          SPAData&   data_end,
    const Int ncols_atm,
    const Int nlevs_atm,
    const Int nswbands,
    const Int nlwbands);

  static void load_and_interp_spa_data(
    const ekat::Comm& comm,
    const std::string&      spa_data_file,
    const SPAInterp&        horiz_weights, 
    const view_1d<Real>&    ps,
          SPAData&          spa_data,
    const Int timelevel,
    const Int nswbands,
    const Int nlwbands);
  

}; // struct Functions

} // namespace spa 
} // namespace scream

#endif // SPA_FUNCTIONS_HPP

// We don't do ETI, since we don't know some of the concrete types.
// E.g., we don't know InputProvider, or ScalarT (although we could
// ETI the "common" cases, where the provider is a view_1d, and
// Scalar=Real or Scalar=Pack<Real,SCREAM_PACK_SIZE>).
# include "spa_functions_impl.hpp"
