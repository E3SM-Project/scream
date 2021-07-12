#include "spa_horizontal_interpolation_impl.hpp"

namespace scream {
namespace spa {

  /*
   * Explicit instantiation for doing update prognostics functions on Reals using the
   * default device.
   */

  template struct SPAFunctions<Real,DefaultDevice>;


} // namespace spa
} // namespace scream
