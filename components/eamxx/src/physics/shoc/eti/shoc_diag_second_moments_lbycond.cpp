#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 *   Explicit instantiation using the default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream

#include "shoc_diag_second_moments_lbycond_impl.hpp"
