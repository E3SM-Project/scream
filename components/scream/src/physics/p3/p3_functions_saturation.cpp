#include "p3_functions_saturation_impl.hpp"
#include "share/scream_types.hpp"

namespace scream {
  namespace p3 {

    /*
     * Explicit instantiation for saturation functions on Reals using the
     * default device.
     */

    template struct Functions<Real,DefaultDevice>;

  } // namespace p3
} // namespace scream
