#ifndef SCREAM_OUTPUT_FILL_POLICIES_HPP
#define SCREAM_OUTPUT_FILL_POLICIES_HPP

#include "share/scream_types.hpp"

namespace scream {
namespace io_policy {

template<typename AvgPolicy>
struct FillDisabled
{
  // Nothing to do when resetting views for Instant average
  template<typename ViewT>
  static void reset_view (ViewT& v) {
    constexpr bool needs_reset_view = AvgPolicy::needs_reset_view;
    if constexpr (needs_reset_view)
      Kokkos::deep_copy (v,AvgPolicy::reset_value);
  }
};

template<typename AvgPolicy>
struct FillEnabled
{
  // Nothing to do when resetting views for Instant average
  template<typename ViewT>
  static void reset_view (ViewT& v) {
    constexpr bool needs_reset_view = AvgPolicy::needs_reset_view;
    if constexpr (needs_reset_view)
      Kokkos::deep_copy (v,DefaultFillValue<Real>::value);
  }
};

} // namespace io_policy
} // namespace scream

#endif // SCREAM_OUTPUT_FILL_POLICIES_HPP
