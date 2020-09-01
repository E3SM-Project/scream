#include "shoc_pblintd_init_pot_impl.hpp"
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

