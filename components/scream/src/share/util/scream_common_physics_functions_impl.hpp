#ifndef SCREAM_COMMON_PHYSICS_IMPL_HPP
#define SCREAM_COMMON_PHYSICS_IMPL_HPP

#include "share/scream_constants.hpp"
#include "share/util/scream_column_ops.hpp"

namespace scream {

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::exner_function(const ScalarT& pressure)
{
  using C = scream::Constants<Real>;

  static constexpr auto p0 = C::P0;
  static constexpr auto rd = C::R_dry_air;
  static constexpr auto inv_cp = C::one / C::cp_dry_air;

  return pow( pressure/p0, rd*inv_cp );
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderP>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::exner_function(const MemberType& team,
                                               const InputProviderP& pressure,
                                               const view_1d<ScalarT>& exner)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,exner.extent(0)),
                       [&] (const int k) {
    exner(k) = exner_function(pressure(k));
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_theta_from_T(const ScalarT& temperature, const ScalarT& pressure)
{
  return temperature/exner_function(pressure);
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderP>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_theta_from_T(const MemberType& team,
                                                       const InputProviderT& temperature,
                                                       const InputProviderP& pressure,
                                                       const view_1d<ScalarT>& theta)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,theta.extent(0)),
                       [&] (const int k) {
    theta(k) = calculate_theta_from_T(temperature(k),pressure(k));
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_T_from_theta(const ScalarT& theta, const ScalarT& pressure)
{
  return theta*exner_function(pressure);
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderP>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_T_from_theta(const MemberType& team,
                                                       const InputProviderT& theta,
                                                       const InputProviderP& pressure,
                                                       const view_1d<ScalarT>& temperature)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,temperature.extent(0)),
                       [&] (const int k) {
    temperature(k) = calculate_T_from_theta(theta(k),pressure(k));
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_temperature_from_virtual_temperature(const ScalarT& T_virtual, const ScalarT& qv)
{
  using C = scream::Constants<Real>;

  // Molecular mass ratio of water vapor and dry air
  static constexpr auto mmr = C::M_water_vapor / C::M_dry_air;
  static constexpr auto one = C::one;

  return T_virtual*((mmr*(one+qv))/(qv+mmr));
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::
calculate_temperature_from_virtual_temperature(const MemberType& team,
                                               const InputProviderT& T_virtual,
                                               const InputProviderQ& qv,
                                               const view_1d<ScalarT>& temperature)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,temperature.extent(0)),
                       [&] (const int k) {
    temperature(k) = calculate_temperature_from_virtual_temperature(T_virtual(k),qv(k));
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::calculate_virtual_temperature(const ScalarT& temperature, const ScalarT& qv)
{
  using C = scream::Constants<Real>;

  // Molecular mass ratio of water vapor and dry air
  static constexpr auto mmr = C::M_water_vapor / C::M_dry_air;
  static constexpr auto one = C::one;

  return temperature*( (qv+mmr)/(mmr*(one+qv)) );
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderQ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::
calculate_virtual_temperature(const MemberType& team,
                              const InputProviderT& temperature,
                              const InputProviderQ& qv,
                              const view_1d<ScalarT>& T_virtual)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,T_virtual.extent(0)),
                       [&] (const int k) {
    T_virtual(k) = calculate_virtual_temperature(temperature(k),qv(k));
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_dse(const ScalarT& temperature, const ScalarT& z, const Real surf_geopotential)
{
  using C = scream::Constants<Real>;

  static constexpr auto cp = C::cp_dry_air;
  static constexpr auto g  = C::gravity;

  return cp*temperature + g*z + surf_geopotential;
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderT, typename InputProviderZ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_dse(const MemberType& team,
                                              const InputProviderT& temperature,
                                              const InputProviderZ& z,
                                              const Real surf_geopotential,
                                              const view_1d<ScalarT>& dse)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,dse.extent(0)),
                       [&] (const int k) {
    dse(k) = calculate_dse(temperature(k),z(k),surf_geopotential);
  });
}

template<typename DeviceT>
template<typename ScalarT>
KOKKOS_INLINE_FUNCTION
ScalarT PhysicsFunctions<DeviceT>::
calculate_dz(const ScalarT& pseudo_density, const ScalarT& p_mid, const ScalarT& T_mid, const ScalarT& qv)
{
  using C = scream::Constants<Real>;

  const ScalarT& T_virtual = calculate_virtual_temperature(T_mid,qv);

  static constexpr auto Rd = C::R_dry_air;
  static constexpr auto g  = C::gravity;
  return (Rd/g)*pseudo_density*T_virtual / p_mid;
}

template<typename DeviceT>
template<typename ScalarT,
         typename InputProviderPD, typename InputProviderP,
         typename InputProviderT,  typename InputProviderQ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_dz(const MemberType& team,
                                             const InputProviderPD& pseudo_density,
                                             const InputProviderP& p_mid,
                                             const InputProviderT& T_mid,
                                             const InputProviderQ& qv,
                                             const view_1d<ScalarT>& dz)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,dz.extent(0)),
                       [&] (const int k) {
    dz(k) = calculate_dz(pseudo_density(k),p_mid(k),T_mid(k),qv(k));
  });
}

template<typename DeviceT>
template<typename ScalarT, typename InputProviderZ>
KOKKOS_INLINE_FUNCTION
void PhysicsFunctions<DeviceT>::calculate_z_int(const MemberType& team,
                                                const int num_levs,
                                                const InputProviderZ& dz,
                                                const Real z_surf,
                                                const view_1d<ScalarT>& z_int)
{
  using column_ops  = ColumnOps<DeviceT,Real>;
  // Note, we set FromTop to false since we are prescribing the *bottom* elevation.
  constexpr bool FromTop = false;
  column_ops::template column_scan<FromTop>(team,num_levs,dz,z_int,z_surf);
}

} // namespace scream

#endif // SCREAM_COMMON_PHYSICS_IMPL_HPP
