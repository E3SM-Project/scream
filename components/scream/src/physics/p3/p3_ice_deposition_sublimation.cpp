#include "p3_ice_deposition_sublimation_impl.hpp"
#include "share/scream_types.hpp"

namespace scream {
namespace p3 {

  /*
   * Explicit instantiation for doing update prognostics functions on Reals using the
   * default device.
   */

  template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
