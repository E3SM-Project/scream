#ifndef SHOC_CONSTANTS_HPP
#define SHOC_CONSTANTS_HPP

//#include "share/scream_types.hpp"

namespace scream {
namespace shoc {

/*
 * Mathematical constants used by shoc.
 *
 * Note that a potential optimization could be to change the type of
 * Scalar constants that have integer values to int.
 */

template <typename Scalar>
struct Constants
{
  static constexpr Scalar mintke = 0.0004; // Minimum TKE [m2/s2]
  static constexpr Scalar maxtke = 50.0;   // Maximum TKE [m2/s2]
};


} // namespace shoc
} // namespace scream

#endif
