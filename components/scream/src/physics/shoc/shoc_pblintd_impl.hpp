#ifndef SHOC_PBLINTD_IMPL_HPP
#define SHOC_PBLINTD_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc pblintd. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::pblintd(const MemberType& team, const Int& nlev, const Int& nlevi, const Int& npbl,
      const uview_1d<const Spack>& z, const uview_1d<const Spack>& zi, 
      const uview_1d<const Spack>& thl, const uview_1d<const Spack>& ql, 
      const uview_1d<const Spack>& q, const uview_1d<const Spack>& u, 
      const uview_1d<const Spack>& v, const uview_1d<const Spack>& cldn,
      const uview_1d<Spack>& rino, const uview_1d<Spack>& thv,
      const Scalar& ustar, const Scalar& obklen, const Scalar& kbfs, 
      Scalar& pblh)
{
  //-----------------------------------------------------------------------
  //
  // Purpose:
  // Diagnose standard PBL variables
  //
  // Method:
  // Diagnosis of PBL depth.
  // The PBL depth follows:
  //    Holtslag, A.A.M., and B.A. Boville, 1993:
  //    Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
  //    Model. J. Clim., vol. 6., p. 1825--1842.
  //
  // Updated by Holtslag and Hack to exclude the surface layer from the
  // definition of the boundary layer Richardson number. Ri is now defined
  // across the outer layer of the pbl (between the top of the surface
  // layer and the pbl top) instead of the full pbl (between the surface and
  // the pbl top). For simiplicity, the surface layer is assumed to be the
  // region below the first model level (otherwise the boundary layer depth
  // determination would require iteration).
  //
  // Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
  // (Use ricr = 0.3 in this formulation)
  //
  // Author: B. Stevens (extracted from pbldiff, August 2000)
  //

  Scalar tlv=0;   // ref. level pot tmp + tmp excess
  bool check;     //  True=>chk if Richardson no.>critcal

  const auto nlev_pack = ekat::npack<Spack>(nlev)-1;
  const int nlev_indx  = (nlev-1)%Spack::n;

   // Compute virtual potential temperature
  shoc_pblintd_init_pot(team, nlev, thl, ql, q, thv);

  pblintd_init(z(nlev_pack)[nlev_indx], check, 
               rino(nlev_pack)[nlev_indx], pblh);

  // PBL height calculation
  pblintd_height(team, nlev, npbl, z, u, v, ustar, thv, 
                 thv(nlev_pack)[nlev_indx], pblh, rino, check);

  // Estimate an effective surface temperature to account for surface
  // fluctuations
  pblintd_surf_temp(nlev, nlevi, npbl, z, ustar, obklen, kbfs,
                    thv, tlv, pblh, check, rino);

  // Improve pblh estimate for unstable conditions using the convective
  // temperature excess as reference temperature:
  pblintd_height(team, nlev, npbl, z, u, v, ustar, thv, tlv, pblh, rino, check);

  // Check PBL height
  pblintd_check_pblh(nlevi, npbl, z, ustar, check, pblh);

  // PBL check over ocean
  shoc_pblintd_cldcheck(zi(nlev_pack)[nlev_indx], cldn(nlev_pack)[nlev_indx], pblh);
}
} // namespace shoc
} // namespace scream

#endif
