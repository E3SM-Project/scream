#ifndef P3_FUNCTIONS_IMPL_HPP
#define P3_FUNCTIONS_IMPL_HPP

#include "p3_functions.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 functions. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */
template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::polysvp1(const Spack& t, const bool ice)
{
  // REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

  // ice
  static Scalar ai[] = {
    6.11147274,     0.503160820,     0.188439774e-1,
    0.420895665e-3, 0.615021634e-5,  0.602588177e-7,
    0.385852041e-9, 0.146898966e-11, 0.252751365e-14};

  // liquid, V1.7
  static Scalar a[] = {
    6.11239921,      0.443987641,     0.142986287e-1,
    0.264847430e-3,  0.302950461e-5,  0.206739458e-7,
    0.640689451e-10,-0.952447341e-13,-0.976195544e-15};

  Spack dt = pack::max(t - 273.16, -80.0);
  Spack result;
  Smask ice_mask = (t < C::Tmelt) && ice;
  Smask liq_mask = !ice_mask;

  // -------------------------------------------

  // Flatau formulation:
  if (ice_mask.any()) {
    Spack ice_result = (ai[0] + dt*(ai[1]+dt*(ai[2]+dt*(ai[3]+dt*(ai[4]+dt*(ai[5]+dt*(ai[6]+dt*(ai[7]+
                                                                                                ai[8]*dt))))))))*100.0;
    result.set(ice_mask, ice_result);
  }
  if (liq_mask.any()) {
    Spack liq_result = (a[0] + dt*(a[1]+dt*(a[2]+dt*(a[3]+dt*(a[4]+dt*(a[5]+dt*(a[6]+dt*(a[7]+a[8]*dt))))))))*100.0;
    result.set(liq_mask, liq_result);
  }

  return result;
}

template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>::qv_sat(const Spack& t_atm, const Spack& p_atm, const bool ice)
{
  Spack e_pres; // saturation vapor pressure [Pa]

  e_pres = polysvp1(t_atm, ice);
  return C::ep_2 * e_pres / pack::max(p_atm-e_pres, 1.e-3);
}

template<typename S, typename D> 
KOKKOS_FUNCTION
void Functions<S,D>
::cloud_water_conservation(const Spack& qc, const Spack& qcnuc,const Scalar dt, 
   Spack& qcaut, Spack& qcacc, Spack &qccol, Spack& qcheti, Spack& qcshd, Spack& qiberg, Spack& qisub, Spack& qidep)
{
  Spack sinks = (qcaut+qcacc+qccol+qcheti+qcshd+qiberg)*dt; // Sinks of cloud water
  Spack sources = qc + (qcnuc)*dt; // Source of cloud water 
  Spack ratio;  

  Smask enforce_conservation  = sinks > sources && sinks >= 1e-20;  // determine if  conservation corrction is necessary 
  Smask nothing_todo = !enforce_conservation; 

  if (enforce_conservation.any()){
    ratio.set(enforce_conservation, sources/sinks);
    qcaut.set(enforce_conservation, qcaut*ratio); 
    qcacc.set(enforce_conservation, qcacc*ratio);
    qccol.set(enforce_conservation, qccol*ratio); 
    qcheti.set(enforce_conservation, qcheti*ratio); 
    qcheti.set(enforce_conservation, qcshd*ratio); 
    qiberg.set(enforce_conservation, qiberg*ratio); 
  }

  if(nothing_todo.any()){
    ratio.set(nothing_todo, 1.0);
  }
  
  enforce_conservation = sources > 1e-20; 
  if (enforce_conservation.any()){
    qidep.set(enforce_conservation, qidep*(1.0-ratio));
    qisub.set(enforce_conservation, qisub*(1.0-ratio));
  }
}

template<typename S, typename D> 
KOKKOS_FUNCTION
void Functions<S,D>
::rain_water_conservation(const Spack& qr, const Spack& qcaut, const Spack& qcacc, const Spack& qimlt, const Spack& qcshd, const Scalar dt,
   Spack& qrevp, Spack& qrcol, Spack& qrheti)
{
  Spack sinks   = (qrevp+qrcol+qrheti)*dt; // Sinks of rain water 
  Spack sources = qr + (qcaut+qcacc+qimlt+qcshd)*dt; // Sources of rain water 
  Spack ratio; 
  
  Smask enforce_conservation  = sinks > sources && sinks >= 1e-20;  // determine if  conservation corrction is necessary 

  if (enforce_conservation.any()){
    sources.set(enforce_conservation, sources/sinks);
    qrevp.set(enforce_conservation, qrevp*ratio); 
    qrcol.set(enforce_conservation, qrcol*ratio);
    qrheti.set(enforce_conservation, qrheti*ratio); 
  }
}

template<typename S, typename D> 
KOKKOS_FUNCTION 
void Functions<S,D>
::ice_water_conservation(const Spack& qitot,const Spack& qidep,const Spack& qinuc,const Spack& qiberg, const Spack &qrcol,const Spack &qccol,const Spack& qrheti,const Spack& qcheti,const Scalar dt, 
   Spack& qisub, Spack& qimlt)  
{
  Spack sinks = (qisub+qimlt)*dt; // Sinks of ice water 
  Spack sources = qitot + (qidep+qinuc+qrcol+qccol+qrheti+qcheti+qiberg)*dt; // Sources of ice water 
  Spack ratio; 

  Smask enforce_conservation  = sinks > sources && sinks >= 1e-20;  // determine if  conservation corrction is necessary 

  if(enforce_conservation.any()){
    sources.set(enforce_conservation, sources/sinks);
    qisub.set(enforce_conservation, qisub*ratio); 
    qimlt.set(enforce_conservation, qimlt*ratio); 
  }
}


} // namespace p3
} // namespace scream

#endif
