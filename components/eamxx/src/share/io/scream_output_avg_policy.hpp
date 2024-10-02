#ifndef SCREAM_OUTPUT_AVG_POLICIES_HPP
#define SCREAM_OUTPUT_AVG_POLICIES_HPP

#include "share/scream_types.hpp"
#include "share/io/scream_io_utils.hpp"

namespace scream {
namespace io_policy {

template<OutputAvgType AvgType, bool UseFill>
struct Combine;

template<bool UseFill>
struct Combine<OutputAvgType::Instant,UseFill>
{
  KOKKOS_INLINE_FUNCTION
  static void apply (const Real& new_val, Real& curr_val)
  {
    if constexpr (UseFill) {
      constexpr Real fill_v = DefaultFillValue<Real>::value;
      if (new_val != fill_v)
        curr_val = new_val;
    } else {
      curr_val = new_val;
    }
  }

  static constexpr bool needs_count = false;
  static constexpr bool needs_reset_view = false;
};

template<bool UseFill>
struct Combine<OutputAvgType::Max,UseFill>
{
  KOKKOS_INLINE_FUNCTION
  static void apply (const Real& new_val, Real& curr_val)
  {
    if constexpr (UseFill) {
      constexpr Real fill_v = DefaultFillValue<Real>::value;
      const bool new_fill  = new_val  == fill_v;
      const bool curr_fill = curr_val == fill_v;
      curr_val = new_fill ? curr_val
                          : (curr_fill ? new_val
                                       : ekat::impl::max(curr_val,new_val));
    } else {
      curr_val = ekat::impl::max(curr_val,new_val);
    }
  }
  static constexpr bool needs_count = false;
  static constexpr bool needs_reset_view = true;
  static constexpr Real reset_value = -std::numeric_limits<Real>::infinity();
};

template<bool UseFill>
struct Combine<OutputAvgType::Min,UseFill>
{
  KOKKOS_INLINE_FUNCTION
  static void apply (const Real& new_val, Real& curr_val)
  {
    if constexpr (UseFill) {
      constexpr Real fill_v = DefaultFillValue<Real>::value;
      const bool new_fill  = new_val  == fill_v;
      const bool curr_fill = curr_val == fill_v;
      curr_val = new_fill ? curr_val
                          : (curr_fill ? new_val
                                       : ekat::impl::min(curr_val,new_val));
    } else {
      curr_val = ekat::impl::min(curr_val,new_val);
    }
  }

  static constexpr bool needs_count = false;
  static constexpr bool needs_reset_view = true;
  static constexpr Real reset_value = std::numeric_limits<Real>::infinity();
};

template<bool UseFill>
struct Combine<OutputAvgType::Average,UseFill>
{
  KOKKOS_INLINE_FUNCTION
  static void apply (const Real& new_val, Real& curr_val)
  {
    if constexpr (UseFill) {
      constexpr Real fill_v = DefaultFillValue<Real>::value;
      const bool new_fill  = new_val  == fill_value;
      const bool curr_fill = curr_val == fill_value;
      curr_val = new_fill ? curr_val
                          : (curr_fill ? new_val
                                       : curr_val+new_val);
    } else {
      curr_val += new_val;
    }
  }

  static constexpr bool needs_count = true;
  static constexpr bool needs_reset_view = true;
  static constexpr Real reset_value = 0;
};

} // namespace io_policy
} // namespace scream

#endif // SCREAM_OUTPUT_AVG_POLICIES_HPP
