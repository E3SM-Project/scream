
#ifndef P3_FUNCTIONS_CHECK_VALUES_IMPL_HPP
#define P3_FUNCTIONS_CHECK_VALUES_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*!------------------------------------------------------------------------------------
  Checks current values of prognotic variables for reasonable values and
  stops and prints values if they are out of specified allowable ranges.
  
  'check_consistency' means include trap for inconsistency in moments;
  otherwise, only trap for Q, T, and negative Qx, etc.  This option is here
  to allow for Q<qsmall.and.N>nsmall or Q>qsmall.and.N<small which can be produced
  at the leading edges due to sedimentation and whose values are accpetable
  since lambda limiters are later imposed after SEDI (so one does not necessarily
  want to trap for inconsistency after sedimentation has been called).
  
  The value 'source_ind' indicates the approximate location in 'p3_main'
  from where 'check_values' was called before it resulted in a trap.
  -----------------------------------------------------------------------------------
*/
template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::check_values(const view_1d<Spack>& qv, const view_1d<Spack>& temp, 
               const Int& kstart, const Int& kend, const MemberType& team, 
               view_1d<Smask> qv_out_bounds, view_1d<Smask> t_out_bounds)
{
  constexpr Scalar T_low  = 173.;
  constexpr Scalar T_high = 323.;
  constexpr Scalar Q_high = 40.e-3;
  constexpr Scalar Q_low  = 0.;

  Int kmin, kmax;
  util::set_min_max(kstart, kend, kmin, kmax, Spack::n);


  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, kmax-kmin+1), [&] (int pk_) {

    const int pk = kmin + pk_;

    const auto t_gt_low_bound   = temp(pk) > T_low;
    const auto t_lt_high_bound   = temp(pk) < T_high;

    const auto qv_gt_low_bound  = qv(pk) > Q_low;
    const auto qv_lt_high_bound = qv(pk) < Q_high;

    t_out_bounds(pk)  = !(t_gt_low_bound && t_lt_high_bound);
    qv_out_bounds(pk) = !(qv_gt_low_bound && qv_lt_high_bound);
  });

}

} // namespace p3
} // namespace scream

#endif
