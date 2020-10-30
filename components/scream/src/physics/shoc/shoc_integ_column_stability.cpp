#include "shoc_integ_column_stability_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing integ_column_stability on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
