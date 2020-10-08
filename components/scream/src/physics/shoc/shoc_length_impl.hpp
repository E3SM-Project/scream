#ifndef SHOC_LENGTH_IMPL_HPP
#define SHOC_LENGTH_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc shoc_length. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_length(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Scalar&                host_dx,
  const Scalar&                host_dy,
  const Scalar&                pblh,
  const uview_1d<const Spack>& tke,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& zi_grid,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& wthv_sec,
  const uview_1d<const Spack>& thv,
  const uview_1d<Spack>&       thv_zi,
  const uview_1d<Spack>&       brunt,
  const uview_1d<Spack>&       shoc_mix)
{
  // Compute the SHOC mixing length scale, which is used to
  // compute the turbulent dissipation in the SGS TKE equation

  // Interpolate virtual potential temperature onto interface grid
  linear_interp(team,zt_grid,zi_grid,thv,thv_zi,nlev,nlevi,0);
  team.team_barrier();

  // Define the brunt vaisalia frequency
  compute_brunt_shoc_length(team,nlev,nlevi,dz_zt,thv,thv_zi,brunt);

  // Find l_inf
  Scalar l_inf = 0;
  compute_l_inf_shoc_length(team,nlev,zt_grid,dz_zt,tke,l_inf);

  // Determine the convective velocity scale of the planetary boundary layer
  Scalar conv_vel = 0;
  compute_conv_vel_shoc_length(team,nlev,pblh,zt_grid,dz_zt,thv,wthv_sec,conv_vel);

  // Compute cubed root of conv_vel and compute eddy turnover timescale.
  Scalar tscale = 0;
  compute_conv_time_shoc_length(pblh,conv_vel,tscale);

  // Compute mixing-length
  compute_shoc_mix_shoc_length(team,nlev,tke,brunt,tscale,zt_grid,l_inf,shoc_mix);

  // Do checks on the length scale.
  check_length_scale_shoc_length(team,nlev,host_dx,host_dy,shoc_mix);
}

} // namespace shoc
} // namespace scream

#endif
