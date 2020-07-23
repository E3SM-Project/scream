#ifndef P3_EVAP_PRECIP_IMPL_HPP
#define P3_EVAP_PRECIP_IMPL_HPP

#include "p3_functions.hpp"
#include "physics_constants.hpp"

namespace scream {
namespace p3 {

/* Evaporate qr in the portion of each cell which has rain but no cloud.
   It is assumed that macrophysics handles condensation/evaporation of qc and
   that there is no condensation of rain. Thus qccon, qrcon and qcevp have
   been removed from the original P3-WRF. Further, sublimation is handled by
   ice_deposition_sublimation since P3 only has a single ice category.
*/ 
template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::evap_precip(
  const Spack& qr_incld, const Spack& qc_incld, const Spack& nr_incld, const Spack& qitot_incld,
  const Spack& lcldm, const Spack& rcldm, const Spack& qvs, const Spack& ab, const Spack& epsr,
  const Spack& qv, Spack& qrevp, Spack& nrevp,
  const Smask& context)
{

  /*Determine temporary cloud fraction, set to zero if cloud water + ice is
    very small.  This will ensure that evap of precip occurs over entire */

  constexpr Scalar QSMALL   = C::QSMALL;


  /*Determine temporary cloud fraction which ramps to 0 as cloud water
    approaches zero. This avoids problems due to lcldm having a positive
    minimum value which would result in evaporation never occurring in a
    small portion of the grid. Ramping rather than thresholding is used
    to preserve convergence.*/
  const auto set_cld_zero = qc_incld < sp(1.e-6) && context;
  Spack cld;
  cld.set(set_cld_zero, lcldm*qc_incld/sp(1.e-6) );
  cld.set(!set_cld_zero && context, lcldm);

  //Only calculate if there is some rain fraction > cloud fraction
  qrevp = 0;
  nrevp = 0;

  const auto do_evap = rcldm > cld && qr_incld >= QSMALL && context;

  if(do_evap.any()){
    //Calculate q for out-of-cloud region
    Spack qclr;
    qclr.set(do_evap,(qv-cld*qvs)/(1-cld));
    //as cld -> 1, qclr can go neg, causing excessive evap. Just limit to zero for now.
    qclr.set(do_evap,pack::max(0,qclr));
    //if qv>qvs (cell-ave supersaturated), qclr can exceedc cell-ave sat and weird things
    //can happen. Not sure it does, but being safe.
    qclr.set(do_evap,pack::min(qclr,qv));

    //Compute evap rate in evaporating region. Note qclr<=qv<=qvs so qclr-qvs
    //should be negative. The minus sign in front makes qrevp positive.  In practice
    //roundoff error *might* make qrevp slightly negative, so adding max(0,...).
    qrevp.set(do_evap, pack::max(0,-epsr * (qclr-qvs)/ab));

    //Turn into cell-average evap rate
    qrevp.set(do_evap,qrevp*(rcldm-cld));

    //and now turn into average over precipitating area
    qrevp.set(do_evap,qrevp/rcldm);

    //reduce drop number proportional to mass lost
    nrevp.set(do_evap,nr_incld*(qrevp/qr_incld));
  }
}


} // namespace p3
} // namespace scream

#endif

