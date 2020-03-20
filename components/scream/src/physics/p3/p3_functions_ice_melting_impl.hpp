#ifndef P3_FUNCTIONS_ICE_MELTING_IMPL_HPP
#define P3_FUNCTIONS_ICE_MELTING_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_melting(const Spack& rho, const Spack& t, const Spack& pres, const Spack& rhofaci,
	      const Spack& f1pr05, const Spack& f1pr14, const Spack& xxlv, const Spack& xlf, 
	      const Spack& dv, const Spack& sc, const Spack& mu, const Spack& kap, 
	      const Spack& qv, const Spack& qitot_incld, const Spack& nitot_incld,
	      Spack& qimlt, Spack& nimlt)
{
  // Notes Left over from WRF Version:
  // need to add back accelerated melting due to collection of ice mass by rain (pracsw1)
  // note 'f1pr' values are normalized, so we need to multiply by N
  // currently enhanced melting from collision is neglected
  // include RH dependence

  const auto Pi = C::Pi;
  const auto QSMALL = C::QSMALL;
  const auto Tmelt = C::Tmelt;
  
  //Find cells above freezing AND which have ice
  const auto has_melt_qi = (qitot_incld >= QSMALL ) && (t > Tmelt);

  if (has_melt_qi.any()){

    //PMC qv_sat from math_impl.hpp seems to match hardcoded formula from F90 I'm swapping in C++ ver.
    //    Note that qsat0 should be with respect to liquid. Confirmed F90 code did this.

    //const auto qsat0 = qv_sat(Spack(Tmelt), pres, false); //last false means NOT saturation w/ respect to ice.
    const auto e0 = polysvp1(Spack(Tmelt), 0);
    const auto qsat0 = 0.622 *e0/(pres-e0);

    qimlt.set(has_melt_qi, ( (f1pr05+f1pr14*pack::cbrt(sc)*pack::sqrt(rhofaci*rho/mu))
			     *((t-Tmelt)*kap-rho*xxlv*dv*(qsat0-qv))
			     *sp(2.0)* Pi /xlf)*nitot_incld );

    //make sure qimlt is always negative
    qimlt = pack::max(qimlt,sp(0.0));
    
    //Reduce ni in proportion to decrease in qi mass. Prev line makes sure it always has the right sign.
    nimlt.set(has_melt_qi, qimlt*(nitot_incld/qitot_incld) );
   
  }
  
}
  
} // namespace p3
} // namespace scream

#endif
