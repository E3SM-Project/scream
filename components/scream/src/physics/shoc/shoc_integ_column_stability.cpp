#include "shoc_integ_column_stability_impl.hpp"
#include "share/scream_types.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
