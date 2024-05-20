#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for using the default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream

#include "shoc_diag_third_shoc_moments_impl.hpp"
