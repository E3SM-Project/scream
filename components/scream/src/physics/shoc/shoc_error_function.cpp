#include "shoc_error_function_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for the default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
