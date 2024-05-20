#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing pblintd_height on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream

#include "shoc_pblintd_height_impl.hpp"
