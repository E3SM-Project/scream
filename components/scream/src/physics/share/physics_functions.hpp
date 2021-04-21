#ifndef PHYSICS_FUNCTIONS_HPP
#define PHYSICS_FUNCTIONS_HPP

#include "physics_constants.hpp"

#include "share/scream_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream {
namespace physics {

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for shared physics. We use the ETI pattern for
 * these functions.
 *
 * Assumptions:
 *  - Kokkos team policies have a vector length of 1
 */

template <typename ScalarT, typename DeviceT>
struct Functions
{

  enum SaturationFcn { Polysvp1 = 0, MurphyKoop = 1};

  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S>
  using BigPack = ekat::Pack<Scalar,SCREAM_PACK_SIZE>;
  template <typename S>
  using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

  using IntSmallPack = SmallPack<Int>;
  using Pack         = BigPack<Scalar>;
  using Spack        = SmallPack<Scalar>;

  template <typename S>
  using Mask = ekat::Mask<Pack::n>;

  template <typename S>
  using SmallMask = ekat::Mask<SmallPack<S>::n>;

  using Smask = SmallMask<Scalar>;

  using KT = KokkosTypes<Device>;

  using C = scream::physics::Constants<Scalar>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;
  template <typename S>
  using uview_2d = typename ekat::template Unmanaged<view_2d<S> >;

  using MemberType = typename KT::MemberType;

  using Workspace = typename ekat::WorkspaceManager<Spack, Device>::Workspace;

  //
  // --------- Functions ---------
  //

  // For all functions, range_mask represents a mask on the Spack variables which
  // determines which elements are active (not masked).

  //  compute saturation vapor pressure
  //  polysvp1 returned in units of pa.
  //  t is input in units of k.
  //  ice refers to saturation with respect to liquid (false) or ice (true)
  KOKKOS_FUNCTION
  static Spack polysvp1(const Spack& t, const bool ice, const Smask& range_mask);

  //  compute saturation vapor pressure using Murphy and Koop(2005) formulation
  //  MurphyKoop_svp returned in units of pa.
  //  t is input in units of k.
  //  ice refers to saturation with respect to liquid (false) or ice (true)
  KOKKOS_FUNCTION
  static Spack MurphyKoop_svp(const Spack& t, const bool ice, const Smask& range_mask);

  // Calls a function to obtain the saturation vapor pressure, and then computes
  // and returns the saturation mixing ratio, with respect to either liquid or ice,
  // depending on value of 'ice'
  KOKKOS_FUNCTION
  static Spack qv_sat(const Spack& t_atm, const Spack& p_atm, const bool ice, const Smask& range_mask, const SaturationFcn func_idx = MurphyKoop);

  //checks temperature for negatives and NaNs
  KOKKOS_FUNCTION
  static void check_temperature(const Spack& t_atm, const char* func_name, const Smask& range_mask);

  // Computes exner function
  // The result is exners formula, and is dimensionless
  // The input is mid-level pressure, and has units of Pa
  KOKKOS_FUNCTION
  static Spack get_exner(const Spack& pmid, const Smask& range_mask);

  // Converts temperature into potential temperature
  // The result is the potential temperature, units in K
  // The inputs are
  //   T_mid is the atmospheric temperature, units in K
  //   p_mid is the atmospheric pressure, units in Pa.  Used with Exners function using `get_exner` defined above.
  KOKKOS_FUNCTION
  static Spack get_potential_temperature(const Spack& T_mid, const Spack& p_mid, const Smask& range_mask);

  template <typename InputProvider>
  KOKKOS_FUNCTION
  static void get_potential_temperature(const int nlev,
                                        const InputProvider& T_mid,
                                        const InputProvider& p_mid,
                                        const view_1d<Scalar>& T_potential);
  // Converts potential temperature into temperature
  // The result is the temperature, units in K
  // The inputs are
  //   Potential Temperature, units in K
  //   p_mid is the atmospheric pressure, units in Pa.  Used with Exners function using `get_exner` defined above.
  KOKKOS_FUNCTION
  static Spack get_potential_temperature_inv(const Spack& th_mid, const Spack& p_mid, const Smask& range_mask);

  template <typename InputProvider>
  KOKKOS_FUNCTION
  static void get_potential_temperature_inv(const int nlev,
                                             const InputProvider& th_mid,
                                             const InputProvider& p_mid,
                                             const view_1d<Scalar>& T_mid);
  // Determine the physical thickness of a vertical layer
  // The result is dz, units in m
  // The inputs are
  //   zi_top is the distance above the surface of the top of the layer.  Units in m
  //   zi_bot is the distance above the surface of the bottom of the layer.  Units in m
  KOKKOS_FUNCTION
  static Spack get_dz(const Spack& zi_top, const Spack& zi_bot, const Smask& range_mask);

  // Compute dry static energy (DSE).
  // The result unit is in J/kg
  // The inputs are
  //   T_mid is the atmospheric temperature. Units in K.
  //   z_mid is the geopotential height above surface at midpoints. Units in m.
  //   surf_geopotential is the surface geopotential height. Units in m.
  KOKKOS_FUNCTION
  static Spack get_dse(const Spack& T_mid, const Spack& z_mid, const Real surf_geopotential, const Smask& range_mask);

  // Compute virtual temperature
  // The result unit is in K
  // The inputs are
  //   T_mid is the atmospheric temperature.  Units in K.
  //   qv    is the water vapor mass mixing ratio.  Units in kg/kg
  KOKKOS_FUNCTION
  static Spack get_virtual_temperature(const Spack& T_mid, const Spack& qv, const Smask& range_mask);

  template <typename InputProvider>
  KOKKOS_FUNCTION
  static void get_virtual_temperature(const int nlev,
                                      const InputProvider& T_mid,
                                      const InputProvider& qv,
                                      const view_1d<Scalar>& T_virtual);
};

} // namespace physics
} // namespace scream

// If a GPU build, make all code available to the translation unit; otherwise,
// ETI is used.
#ifdef KOKKOS_ENABLE_CUDA
# include "physics_saturation_impl.hpp"
# include "physics_universal_impl.hpp"
#endif

#endif // PHYSICS_FUNCTIONS_HPP
