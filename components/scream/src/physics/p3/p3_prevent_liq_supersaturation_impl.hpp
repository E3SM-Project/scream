#ifndef P3_PREVENT_LIQ_SUPERSATURATION_IMPL_HPP
#define P3_PREVENT_LIQ_SUPERSATURATION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_functions.hpp" // also for ETI not on GPUs
#include "physics/share/physics_saturation_impl.hpp"

namespace scream {
namespace p3 {

//Limit sublimation and evaporation to ensure qv doesn't become supersaturated.
//This becomes a bit subtle because of the difference between condensational
//versus sublimational heating in the psychrometric correction. 
  
template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::prevent_liq_supersaturation(const Spack& pres, const Spack& t_atm, const Spack& qv, const Spack& latent_heat_vapor, const Spack& latent_heat_sublim, const Scalar& dt, const Spack& qv2qi_vapdep_tend, const Spack& qinuc, Spack& qi2qv_sublim_tend, Spack& qr2qv_evap_tend, const Smask& range_mask,
  const Smask& context)
// Note: range_mask handles padded data in packed arrays and context handles legit cells that we
// know don't have any condensate worth acting on
{
  using physics = scream::physics::Functions<Scalar, Device>;
  
  constexpr Scalar inv_cp       = C::INV_CP;
  constexpr Scalar rv           = C::RV;

  Spack qv_sinks, qv_sources, qv_endstep, T_endstep, A, frac;
  
  qv_sinks.set(context, qv2qi_vapdep_tend + qinuc);
  qv_sources.set(context, qi2qv_sublim_tend + qr2qv_evap_tend);

  //printf("qv=%f\n",qv[0]);
  //printf("qv_sources=%f\n",qv_sources[0]);
  //printf("t_atm=%f\n",t_atm[0]);

  //printf("dt=%f\n",dt);
  
  //Actual qv and T after microphys step
  qv_endstep.set(context,qv - qv_sinks*dt + qv_sources*dt);
  T_endstep.set(context,t_atm + ( (qv_sinks-qi2qv_sublim_tend)*latent_heat_sublim*inv_cp
				  - qr2qv_evap_tend*latent_heat_vapor*inv_cp )*dt);

  //printf("qv_endstep=%f\n",qv_endstep[0]);
  //printf("T_endstep=%f\n",T_endstep[0]);
  
  //qv we would have at end of step if we were saturated with respect to liquid
  const auto qsl = physics::qv_sat(T_endstep,pres,false,range_mask); //"false" means NOT sat w/ respect to ice

  //printf("qsl=%f\n",qsl[0]);
  
  //The balance we seek is:
  // qv-qv_sinks*dt+qv_sources*frac*dt=qsl+dqsl_dT*(T correction due to conservation)
  // where the T correction for conservation is:
  // dt*[latent_heat_sublim/cp*(qi2qv_sublim_tend-frac*qi2qv_sublim_tend)
  //     +latent_heat_vapor/cp*(qr2qv_evap_tend  -frac*qr2qv_evap_tend)]
  // =(1-frac)*dt/cp*(latent_heat_sublim*qi2qv_sublim_tend + latent_heat_vap*qr2qv_evap_tend).
  // Note T correction is positive because frac *reduces* evaporative cooling. Note as well that
  // dqsl_dt comes from linearization of qsl around the end-of-step T computed before temperature
  // correction. dqsl_dt should be computed with respect to *liquid* even though frac also adjusts
  // sublimation because we want to be saturated with respect to liquid at the end of the step.
  // dqsl_dt=Latent_heat_vapor*qsl/rv*T^2 following Clausius Clapeyron. Combining and solving for
  // frac yields:
   
  A.set(context,latent_heat_vapor*qsl*dt*inv_cp/(rv*T_endstep*T_endstep)
	* (latent_heat_sublim*qi2qv_sublim_tend + latent_heat_vapor*qr2qv_evap_tend) );

  
  frac.set(context, (qsl-qv+qv_sinks*dt + A)/(qv_sources*dt + A) );

  //The only way frac<0 is if qv-qv_sinks*dt is already greater than qsl. In this case
  //the best we can do is zero out qv_sources.
  frac = max(0,frac);

  //The only way frac>1 is if qv-qv_sinks*dt+qv_sources*dt < qsl, in which case we shouldn't 
  //limit anyways so set frac to 1:
  frac = min(1,frac);

  //printf("A=%f\n",A[0]);
  //printf("frac=%f\n",frac[0]);
  
  qi2qv_sublim_tend.set(context, frac*qi2qv_sublim_tend );
  qr2qv_evap_tend.set(context, frac*qr2qv_evap_tend);

} //end of prevent_liq_supersaturation

} // namespace p3
} // namespace scream

#endif
