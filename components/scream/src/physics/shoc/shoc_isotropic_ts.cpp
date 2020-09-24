#include "shoc_isotropic_ts_impl.hpp"
#include "share/scream_types.hpp"

namespace scream {
  namespace shoc {

    /*
     * Explicit instantiation using the default device.
     */

    template struct Functions<Real,DefaultDevice>;

  } // namespace shoc
} // namespace scream
