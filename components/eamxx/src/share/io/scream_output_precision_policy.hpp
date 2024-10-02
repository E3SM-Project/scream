#ifndef SCREAM_OUTPUT_PRECISION_POLICIES_HPP
#define SCREAM_OUTPUT_PRECISION_POLICIES_HPP

#include "share/scream_types.hpp"

namespace scream {
namespace io_policy {

struct SinglePrecision
{
  static constexpr const char* name = "float";
  using type = float;
}

struct SinglePrecision
{
  static constexpr const char* name = "float";
  using type = double;
}

} // namespace io_policy
} // namespace scream

#endif // SCREAM_OUTPUT_PRECISION_POLICIES_HPP

