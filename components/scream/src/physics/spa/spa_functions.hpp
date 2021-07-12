#ifndef SPA_FUNCTIONS_HPP
#define SPA_FUNCTIONS_HPP

#include "share/scream_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

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

  using Mask = ekat::Mask<BigPack<Scalar>::n>;
  using Smask = ekat::Mask<SmallPack<Scalar>::n>;

  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using view_3d = typename KT::template view_3d<S>;
  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

  struct PhysicalState {
    PhysicalState() = default;
    // Surface and Reference Pressure:
  };

  struct MonthlyGHG {
    MonthlyGHG() = default;
    // Aerosol shortwave assymmetry parameter
    view_3d<Spack> AER_G_SW_0;
    view_3d<Spack> AER_G_SW_1;
    view_3d<Spack> AER_G_SW_2;
    view_3d<Spack> AER_G_SW_3;
    view_3d<Spack> AER_G_SW_4;
    view_3d<Spack> AER_G_SW_5;
    view_3d<Spack> AER_G_SW_6;
    view_3d<Spack> AER_G_SW_7;
    view_3d<Spack> AER_G_SW_8;
    view_3d<Spack> AER_G_SW_9;
    view_3d<Spack> AER_G_SW_10;
    view_3d<Spack> AER_G_SW_11;
    view_3d<Spack> AER_G_SW_12;
    view_3d<Spack> AER_G_SW_13;
    // Aerosol shortwave single scattering albedo
    view_3d<Spack> AER_SSA_SW_0;
    view_3d<Spack> AER_SSA_SW_1;
    view_3d<Spack> AER_SSA_SW_2;
    view_3d<Spack> AER_SSA_SW_3;
    view_3d<Spack> AER_SSA_SW_4;
    view_3d<Spack> AER_SSA_SW_5;
    view_3d<Spack> AER_SSA_SW_6;
    view_3d<Spack> AER_SSA_SW_7;
    view_3d<Spack> AER_SSA_SW_8;
    view_3d<Spack> AER_SSA_SW_9;
    view_3d<Spack> AER_SSA_SW_10;
    view_3d<Spack> AER_SSA_SW_11;
    view_3d<Spack> AER_SSA_SW_12;
    view_3d<Spack> AER_SSA_SW_13;
    // Aerosol longwave absorption optical depth
    view_3d<Spack> AER_TAU_LW_0;
    view_3d<Spack> AER_TAU_LW_1;
    view_3d<Spack> AER_TAU_LW_2;
    view_3d<Spack> AER_TAU_LW_3;
    view_3d<Spack> AER_TAU_LW_4;
    view_3d<Spack> AER_TAU_LW_5;
    view_3d<Spack> AER_TAU_LW_6;
    view_3d<Spack> AER_TAU_LW_7;
    view_3d<Spack> AER_TAU_LW_8;
    view_3d<Spack> AER_TAU_LW_9;
    view_3d<Spack> AER_TAU_LW_10;
    view_3d<Spack> AER_TAU_LW_11;
    view_3d<Spack> AER_TAU_LW_12;
    view_3d<Spack> AER_TAU_LW_13;
    view_3d<Spack> AER_TAU_LW_14;
    view_3d<Spack> AER_TAU_LW_15;
    // Aerosol shortwave extinction optical depth
    view_3d<Spack> AER_TAU_SW_0;
    view_3d<Spack> AER_TAU_SW_1;
    view_3d<Spack> AER_TAU_SW_2;
    view_3d<Spack> AER_TAU_SW_3;
    view_3d<Spack> AER_TAU_SW_4;
    view_3d<Spack> AER_TAU_SW_5;
    view_3d<Spack> AER_TAU_SW_6;
    view_3d<Spack> AER_TAU_SW_7;
    view_3d<Spack> AER_TAU_SW_8;
    view_3d<Spack> AER_TAU_SW_9;
    view_3d<Spack> AER_TAU_SW_10;
    view_3d<Spack> AER_TAU_SW_11;
    view_3d<Spack> AER_TAU_SW_12;
    view_3d<Spack> AER_TAU_SW_13;
  };

  struct MonthlyCCN {
    MonthlyCCN() = default;
    // CCN concentration at S=0.1%
    view_3d<Spack> CCN;
  };

  struct PrescribedAero {
    PrescribedAero() = default;
    // CCN concentration at S=0.1%
    view_2d<Spack> CCN;
    // Aerosol shortwave assymmetry parameter
    view_2d<Spack> AER_G_SW_0;
    view_2d<Spack> AER_G_SW_1;
    view_2d<Spack> AER_G_SW_2;
    view_2d<Spack> AER_G_SW_3;
    view_2d<Spack> AER_G_SW_4;
    view_2d<Spack> AER_G_SW_5;
    view_2d<Spack> AER_G_SW_6;
    view_2d<Spack> AER_G_SW_7;
    view_2d<Spack> AER_G_SW_8;
    view_2d<Spack> AER_G_SW_9;
    view_2d<Spack> AER_G_SW_10;
    view_2d<Spack> AER_G_SW_11;
    view_2d<Spack> AER_G_SW_12;
    view_2d<Spack> AER_G_SW_13;
    // Aerosol shortwave single scattering albedo
    view_2d<Spack> AER_SSA_SW_0;
    view_2d<Spack> AER_SSA_SW_1;
    view_2d<Spack> AER_SSA_SW_2;
    view_2d<Spack> AER_SSA_SW_3;
    view_2d<Spack> AER_SSA_SW_4;
    view_2d<Spack> AER_SSA_SW_5;
    view_2d<Spack> AER_SSA_SW_6;
    view_2d<Spack> AER_SSA_SW_7;
    view_2d<Spack> AER_SSA_SW_8;
    view_2d<Spack> AER_SSA_SW_9;
    view_2d<Spack> AER_SSA_SW_10;
    view_2d<Spack> AER_SSA_SW_11;
    view_2d<Spack> AER_SSA_SW_12;
    view_2d<Spack> AER_SSA_SW_13;
    // Aerosol longwave absorption optical depth
    view_2d<Spack> AER_TAU_LW_0;
    view_2d<Spack> AER_TAU_LW_1;
    view_2d<Spack> AER_TAU_LW_2;
    view_2d<Spack> AER_TAU_LW_3;
    view_2d<Spack> AER_TAU_LW_4;
    view_2d<Spack> AER_TAU_LW_5;
    view_2d<Spack> AER_TAU_LW_6;
    view_2d<Spack> AER_TAU_LW_7;
    view_2d<Spack> AER_TAU_LW_8;
    view_2d<Spack> AER_TAU_LW_9;
    view_2d<Spack> AER_TAU_LW_10;
    view_2d<Spack> AER_TAU_LW_11;
    view_2d<Spack> AER_TAU_LW_12;
    view_2d<Spack> AER_TAU_LW_13;
    view_2d<Spack> AER_TAU_LW_14;
    view_2d<Spack> AER_TAU_LW_15;
    // Aerosol shortwave extinction optical depth
    view_2d<Spack> AER_TAU_SW_0;
    view_2d<Spack> AER_TAU_SW_1;
    view_2d<Spack> AER_TAU_SW_2;
    view_2d<Spack> AER_TAU_SW_3;
    view_2d<Spack> AER_TAU_SW_4;
    view_2d<Spack> AER_TAU_SW_5;
    view_2d<Spack> AER_TAU_SW_6;
    view_2d<Spack> AER_TAU_SW_7;
    view_2d<Spack> AER_TAU_SW_8;
    view_2d<Spack> AER_TAU_SW_9;
    view_2d<Spack> AER_TAU_SW_10;
    view_2d<Spack> AER_TAU_SW_11;
    view_2d<Spack> AER_TAU_SW_12;
    view_2d<Spack> AER_TAU_SW_13;
  };

  struct InterpolationData {
    InterpolationData() = default;
    // Needed for temporal interpolation
    Real                 current_time_in_days;
    Real                 length_of_month_in_days;
    // Needed for vertical interpolation
    view_2d<const Spack> p_mid;
  };

  static void init();

  static void main(
    const Int nj, 
    const Int nk,
    const Int c_month,
    const InterpolationData InterpData,
    const MonthlyGHG& GHG,
    const MonthlyCCN& CCN,
    const PrescribedAero& prescribed_aero);

  KOKKOS_FUNCTION
  static void spatial_interpolation( 
    const MemberType& team,
    const Int& nk,
    const uview_1d<const Spack>& c_month_data,
    const uview_1d<const Spack>& n_month_data,
    const uview_1d<Spack>& c_month_interp_data,
    const uview_1d<Spack>& n_month_interp_data
  );

  KOKKOS_FUNCTION
  static void horizontal_interpolation(const std::string& filename);

}; // struct SPAFunctions

} // namespace spa 
} // namespace scream

#endif // SPA_FUNCTIONS_HPP
