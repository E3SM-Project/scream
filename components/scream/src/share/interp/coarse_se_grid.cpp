#include <share/interp/coarse_se_grid.hpp>

namespace scream {
namespace interpolators {

CoarseSEGrid::CoarseSEGrid(int ne_, int np_):
  ne(ne_), np(np_) {
  EKAT_ASSERT(ne_ > 0);
  EKAT_ASSERT(np_ >= 2);
  EKAT_ASSERT((np_ % 2) == 0);
  EKAT_ASSERT(np_ <= this->np_max);

}

} // namespace interpolators
} // namespace scream


