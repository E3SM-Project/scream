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
::calc_bulkRhoRime(const Spack& qi_tot, Spack& qi_rim, Spack& bi_rim, Spack& rho_rime){

auto rim_exists = bi_rim > 1e-15; 
auto rim_doesnot_exist = !rim_exists; 
if (rim_exists.any()){
  bi_rim.set(rim_exists, qi_rim/bi_rim); 
  // impose limits on rho_rime;  adjust bi_rim if needed
  auto small_rho_rime = rho_rime < C::rho_rimeMin;  
  auto not_small_rho_rime = !small_rho_rime; 
  if(small_rho_rime.any()){
    rho_rime.set(small_rho_rime, C::rho_rimeMin);
    bi_rim.set(small_rho_rime, qi_rim/rho_rime);
  }

  if(not_small_rho_rime.any()){
    rho_rime.set(not_small_rho_rime, C::rho_rimeMax);
    bi_rim.set(not_small_rho_rime, qi_rim/rho_rime);
  }
}

if(rim_doesnot_exist.any()){
  qi_rim.set(rim_doesnot_exist, 0.0); 
  bi_rim.set(rim_doesnot_exist, 0.0);
  rho_rime.set(rim_doesnot_exist, 0.0); 
}

// set upper constraint qi_rim <= qi_tot
auto upper_limit_qi_rim = qi_rim > qi_tot && rho_rime > 0.0; 
if(upper_limit_qi_rim.any()){
 qi_rim.set(upper_limit_qi_rim, qi_tot); 
 bi_rim.set(upper_limit_qi_rim, qi_rim/rho_rime);  
}

//impose consistency 
auto enforce_consistency = qi_rim < C::QSMALL; 
if(enforce_consistency.any()){
  qi_rim.set(enforce_consistency, 0.0);
  bi_rim.set(enforce_consistency, 0.0); 
}

}

template<typename S, typename D> 
KOKKOS_FUNCTION
void Functions<S,D>
::impose_max_total_Ni(Spack& nitot_local, Spack& inv_rho_local){

  auto impose_maximum = nitot_local >= 1e-20; 
  if(impose_maximum.any()){
    auto nitot_imposed_max = C::max_total_Ni * inv_rho_local/nitot_local; 
    nitot_local.set(impose_maximum, nitot_local * min(nitot_imposed_max, 1.0)); 
  }

}

template<typename S, typename D> 
KOKKOS_FUNCTION
void Functions<S,D>
// Khroutdinov and Kogan (2000). The Fortran version supports other options by setting iparam. 
// Here in the C++ code we are using iparam=3 option. 
::cloud_water_autoconversion(const Spack& rho, const Spack& qc_incld,
    const Spack& nc_incld, Spack& qcaut, 
    Spack& ncautc, Spack& ncautr)
{
  auto qc_not_small = qc_incld >= 1e-8; 
  if(qc_not_small.any()){
    qcaut.set(qc_not_small, 1350.0*pow(qc_incld,2.47)*pow(nc_incld*1.e-6*rho,-1.79)); 
    ncautr.set(qc_not_small, qcaut * C::CONS3);
    ncautc.set(qc_not_small, qcaut*nc_incld/qc_incld);  
  }
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
    ratio.set(nothing_todo, 1.0); // If not limiting sinks on qc then most likely did not run out of qc
  }
  
  
  //PMC: ratio is also frac of step w/ liq. thus we apply qiberg for
  //"ratio" of timestep and vapor deposition and sublimation  for the 
  //remaining frac of the timestep.  Only limit if there will be cloud
  // water to begin with.
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

template<typename S, typename D> 
KOKKOS_FUNCTION 
void Functions<S,D>
::  back_to_cell_average(const Spack& lcldm,const Spack& rcldm, const Spack& icldm,    
   Spack& qcacc, Spack& qrevp, Spack& qcaut, Spack& ncacc, Spack& ncslf, Spack& ncautc, Spack& nrslf, 
   Spack& nrevp, Spack& ncautr, Spack& qcnuc,Spack& ncnuc, Spack& qisub, Spack& nrshdr, Spack& qcheti, 
   Spack& qrcol, Spack& qcshd, Spack& qimlt, Spack& qccol, Spack& qrheti, Spack& nimlt, Spack& nccol,
   Spack& ncshdc, Spack& ncheti, Spack& nrcol, Spack& nislf, Spack& qidep, Spack& nrheti, Spack& nisub,
   Spack& qinuc, Spack& ninuc, Spack& qiberg){

     Spack ir_cldm(min(icldm, rcldm)); // Intersection of ICE and RAIN cloud
     Spack il_cldm(min(icldm, lcldm)); // Intersection of ICE and LIQUID cloud
     Spack lr_cldm(min(lcldm, rcldm)); // Intersection of LIQUID and RAIN cloud
     
     //Some process rates take place within the intersection of liquid, rain and ice cloud fractions.
     //We calculate the intersection as the minimum between combinations of cloud fractions and use 
     //these values to map back to cell-average quantities where applicable.

     //map warm-phase process rates to cell-avg
     qcacc   = qcacc*lr_cldm;     // Accretion of liquid to rain
     qrevp   = qrevp*rcldm;       // Evaporation of rain
     qcaut   = qcaut*lcldm;       // Autoconversion of liquid
     ncacc   = ncacc*lr_cldm;     // Number change due to accretion
     ncslf   = ncslf*lcldm;       // Self collection occurs locally in liq. cloud
     ncautc  = ncautc*lcldm;      // Impact of autoconversion on number
     nrslf   = nrslf*rcldm;       // Self collection occurs locally in rain cloud
     nrevp   = nrevp*rcldm;       // Change in rain number due to evaporation 
     ncautr  = ncautr*lr_cldm;    // Autoconversion of rain drops within rain/liq cloud
     
     // AaronDonahue: These variables are related to aerosol activation and their usage will be changed in a later PR.
     qcnuc   = qcnuc*lcldm;       // Impact on liq. from nucleation
     ncnuc   = ncnuc*lcldm;       // Number change due to aerosol activation
     
     // map ice-phase  process rates to cell-avg
     qisub   = qisub*icldm;       // Sublimation of ice in ice cloud
     nrshdr  = nrshdr*il_cldm;    //  Rain # increase due to shedding from rain-ice collisions, occurs when ice and liquid interact
     qcheti  = qcheti*il_cldm;    // Immersion freezing of cloud drops
     qrcol   = qrcol*ir_cldm;     // Collection of rain mass by ice
     qcshd   = qcshd*il_cldm;     // Rain mass growth due to shedding of fain drops after collisions with ice, occurs when ice and liquid interact
     qimlt   = qimlt*icldm;       // Melting of ice
     qccol   = qccol*il_cldm;     // Collection of water by ice
     qrheti  = qrheti*rcldm;      // Immersion freezing of rain
     nimlt   = nimlt*icldm;       // Change in number due to melting
     nccol   = nccol*il_cldm;     // Cloud # change due to collection of cld water by ice
     ncshdc  = ncshdc*il_cldm;    // Number change due to shedding, occurs when ice and liquid interact
     ncheti  = ncheti*lcldm;      // Number change associated with freexzing of cld drops
     nrcol   = nrcol*ir_cldm;     // Rain number change due to collection from ice
     nislf   = nislf*icldm;       // Ice self collection
     qidep   = qidep*icldm;       // Vapor deposition to ice phase
     nrheti  = nrheti*rcldm;      // Change in number due to immersion freezing of rain
     nisub   = nisub*icldm;       // Number change due to sublimation of ice
     
     qiberg  = qiberg*il_cldm;    // Bergeron process
     
     // AaronDonahue: These variables are related to aerosol activation and their usage will be changed in a later PR.
     qinuc   = qinuc;             // Deposition and condensation-freezing nucleation, already cell-averaged
     ninuc   = ninuc;             // Number change due to deposition and condensation-freezing, already cell-averaged 
   }

} // namespace p3
} // namespace scream

#endif
