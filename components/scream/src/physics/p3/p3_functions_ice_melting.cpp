#include "p3_functions_ice_melting_impl.hpp"
#include "share/scream_types.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing p3 ice melting on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
