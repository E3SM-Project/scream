#include "shoc_diag_second_moments_ubycond_impl.hpp"
#include "ekat/scream_types.hpp"

namespace scream {
namespace shoc {

/*
 *  * Explicit instantiation for doing ice melting on Reals using the
 *   * default device.
 *    */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream

