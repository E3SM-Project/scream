#ifndef SCREAM_FIELD_MONOTONICITY_CHECK_HPP
#define SCREAM_FIELD_MONOTONICITY_CHECK_HPP

#include "share/field/field.hpp"

namespace scream
{

// This field property check returns true if all data in a given field are
// monotonically increasing or decreasing, and false if not. Currently, no
// repairs to nonmonotonic fields can be performed. There are ways of repairing
// non-monotonic fields given certain assumptions, but we do not make such
// assumptions here.
template<typename RealType>
class FieldMonotonicityCheck: public FieldPropertyCheck<RealType> {
public:
  using non_const_RT = typename FieldPropertyCheck<RealType>::non_const_RT;
  using const_RT     = typename FieldPropertyCheck<RealType>::const_RT;

  // Default constructor.
  FieldMonotonicityCheck () {}

  // Overrides.

  // The name of the field check
  std::string name () const { return "Monotonicity Field Check"; }

  bool check(const Field<const_RT>& field) const override {
    auto view = field.get_view();
    non_const_RT sign;
    Kokkos::parallel_reduce(view.extent(0), KOKKOS_LAMBDA(Int i, non_const_RT& s) {
      if ((i > 0) && (i < view.extent_int(0)-1)) {
        auto diff1 = view(i) - view(i-1);
        auto diff2 = view(i+1) - view(i);
        s *= (diff1 * diff2 > 0) ? (diff1 * diff2) : 0;
      } else {
        s *= 1;
      }
    }, Kokkos::Prod<non_const_RT>(sign));
    return (sign > 0);
  }

  bool can_repair() const override {
    return false;
  }

  void repair(Field<RealType>& /* field */) const override {
    EKAT_REQUIRE_MSG(false, "Cannot repair a non-monotonic field!");
  }
};

} // namespace scream

#endif // SCREAM_FIELD_MONOTONICITY_CHECK_HPP
