#ifndef SHOC_DIAG_SECOND_SHOC_MOMENTS_IMPL_HPP
#define SHOC_DIAG_SECOND_SHOC_MOMENTS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc diag_second_shoc_moments. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::diag_second_shoc_moments(const MemberType& team, const Int& nlev, const Int& nlevi, 
       const uview_1d<const Spack>& thetal, const uview_1d<const Spack>& qw, const uview_1d<const Spack>& u_wind, 
       const uview_1d<const Spack>& v_wind, const uview_1d<const Spack>& tke, const uview_1d<const Spack>& isotropy,
       const uview_1d<const Spack>& tkh, const uview_1d<const Spack>& tk, const uview_1d<const Spack>& dz_zi, 
       const uview_1d<const Spack>& zt_grid, const uview_1d<const Spack>& zi_grid, const uview_1d<const Spack>& shoc_mix, 
       const Scalar& wthl_sfc, const Scalar& wqw_sfc, const Scalar& uw_sfc, const Scalar& vw_sfc, Scalar& ustar2, Scalar& wstar, 
       const uview_1d<Spack>& isotropy_zi, const uview_1d<Spack>& tkh_zi, const uview_1d<Spack>& tk_zi, const uview_1d<Spack>& thl_sec, 
       const uview_1d<Spack>& qw_sec, const uview_1d<Spack>& wthl_sec, const uview_1d<Spack>& wqw_sec, const uview_1d<Spack>& qwthl_sec, 
       const uview_1d<Spack>& uw_sec, const uview_1d<Spack>& vw_sec, const uview_1d<Spack>& wtke_sec, const uview_1d<Spack>& w_sec)
{
  // This is the main routine to compute the second
  // order moments in SHOC.

  // Calculate surface properties needed for lower
  //  boundary conditions
  shoc_diag_second_moments_srf(wthl_sfc, uw_sfc, vw_sfc, ustar2, wstar);
  team.team_barrier();

  // Diagnose the second order moments flux, for the lower boundary
  const auto nlevi_pack = ekat::npack<Spack>(nlevi)-1;
  const int nlevi_indx  = (nlevi-1)%Spack::n;
  Scalar& wthl_s         = wthl_sec(nlevi_pack)[nlevi_indx];
  Scalar& wqw_s          = wqw_sec(nlevi_pack)[nlevi_indx];
  Scalar& uw_s           = uw_sec(nlevi_pack)[nlevi_indx];
  Scalar& vw_s           = vw_sec(nlevi_pack)[nlevi_indx];
  Scalar& wtke_s         = wtke_sec(nlevi_pack)[nlevi_indx];
  Scalar& thl_s          = thl_sec(nlevi_pack)[nlevi_indx];
  Scalar& qw_s           = qw_sec(nlevi_pack)[nlevi_indx];
  Scalar& qwthl_s        = qwthl_sec(nlevi_pack)[nlevi_indx];
  shoc_diag_second_moments_lbycond(wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, ustar2, wstar,
                                   wthl_s, wqw_s, uw_s, vw_s, wtke_s, thl_s, qw_s, qwthl_s);
  team.team_barrier();

  // Diagnose the second order moments, for points away from boundaries.  this is
  //  the main computation for the second moments
  diag_second_moments(team, nlev, nlevi,
                     thetal, qw, u_wind,v_wind, tke, isotropy,tkh, tk, dz_zi, zt_grid, zi_grid, shoc_mix,
                     isotropy_zi, tkh_zi, tk_zi, thl_sec, qw_sec, wthl_sec, wqw_sec,
                     qwthl_sec, uw_sec, vw_sec, wtke_sec, w_sec);
  team.team_barrier();

  // Diagnose the second order moments, calculate the upper boundary conditions
  Scalar& thl_s0   = thl_sec(0)[0];
  Scalar& qw_s0    = qw_sec(0)[0];
  Scalar& wthl_s0  = wthl_sec(0)[0];
  Scalar& wqw_s0   = wqw_sec(0)[0];
  Scalar& qwthl_s0 = qwthl_sec(0)[0];
  Scalar& uw_s0    = uw_sec(0)[0];
  Scalar& vw_s0    = vw_sec(0)[0];
  Scalar& wtke_s0  = wtke_sec(0)[0];
  shoc_diag_second_moments_ubycond(thl_s0, qw_s0, wthl_s0, wqw_s0, qwthl_s0, uw_s0, vw_s0, wtke_s0);
}

} // namespace shoc
} // namespace scream

#endif
