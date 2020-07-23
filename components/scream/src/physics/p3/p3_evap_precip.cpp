#include "p3_evap_precip_impl.hpp"
#include "ekat/scream_types.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing p3 rain evap on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
