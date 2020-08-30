
#ifndef SHOC_DIAG_SECOND_MOMENTS_LBYCOND_IMPL_HPP
#define SHOC_DIAG_SECOND_MOMENTS_LBYCOND_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "physics_functions.hpp" // also for ETI not on GPUs

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_diag_second_moments_lbycond(
    const Scalar& wthl_sfc, const Scalar& wqw_sfc, const Scalar& uw_sfc, const Scalar& vw_sfc,
    const Scalar& ustar2, const Scalar& wstar,
    Scalar& wthl_sec, Scalar& wqw_sec, Scalar& uw_sec, Scalar& vw_sec,
    Scalar& wtke_sec, Scalar& thl_sec, Scalar& qw_sec, Scalar& qwthl_sec)
{

//printf("c++-1 %16.9f, %16.9f, %16.9f, %16.9f, %16.9f, %16.9f, %16.9f\n",wthl_sfc(0)[0],wqw_sfc(0)[0],ustar2(0)[0],wstar(0)[0],thl_sec(0)[0],uw_sfc(0)[0],vw_sfc(0)[0]);

  // Purpose of this subroutine is to diagnose the lower
  // boundary condition for the second order moments needed
  // for the SHOC parameterization.
  // The thermodymnamic, tracer, and momentum fluxes are set
  // to the surface fluxes for the host model, while the
  // thermodynamic variances and covariances are computed
  // according to that of Andre et al. 1978.

  const auto one_p_eight = 1.8;
  const auto ufmin = 0.01;
  const auto zero_p_two   = 0.2;
  const auto zero_p_three = 0.3;
  const auto zero_p_four  = 0.4;
  const auto two   = 2.0;
  const auto three = 3.0;


  auto uf = sqrt(ustar2+zero_p_three*wstar*wstar);
  uf = max(ufmin,uf);

  // Diagnose thermodynamics variances and covariances
  thl_sec   = zero_p_four*one_p_eight*pow(wthl_sfc/uf, two);
  qw_sec    = zero_p_four*one_p_eight*pow(wqw_sfc/uf, two); 
  qwthl_sec = zero_p_two*one_p_eight*(wthl_sfc/uf)*(wqw_sfc/uf);

  // Vertical fluxes of heat and moisture, simply
  // use the surface fluxes given by host model
  wthl_sec = wthl_sfc;
  wqw_sec  = wqw_sfc;
  uw_sec   = uw_sfc;
  vw_sec   = vw_sfc;
  wtke_sec = pow(max(sqrt(ustar2),ufmin), three);

}


} // namespace shoc
} // namespace scream

#endif
