#include "shoc_eddy_diffusivities_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing eddy_diffusivities on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
