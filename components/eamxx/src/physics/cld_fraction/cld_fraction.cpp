#include "physics/cld_fraction/cld_fraction_main_impl.hpp"
#include "share/core/scream_types.hpp"

namespace scream {
namespace cld_fraction {

/*
 * Explicit instantiation for doing main on Reals using the
 * default device.
 */

template struct CldFractionFunctions<Real,DefaultDevice>;

} // namespace cld_fraction
} // namespace scream
