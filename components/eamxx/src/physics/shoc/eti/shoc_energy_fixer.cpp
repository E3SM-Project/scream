#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation using the default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream

#include "shoc_energy_fixer_impl.hpp"
