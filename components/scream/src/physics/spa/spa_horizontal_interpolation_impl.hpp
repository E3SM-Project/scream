#ifndef SPA_HORIZONTAL_INTERPOLATION_IMPL_HPP
#define SPA_HORIZONTAL_INTERPOLATION_IMPL_HPP

#include "physics/spa/spa_functions.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/util/ekat_string_utils.hpp"

namespace scream {
namespace spa {

/*
 *  TODO: Add comments here on what the horizontal interpolation is doing.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void SPAFunctions<S,D>
::horizontal_interpolation(const std::string& remap_filename)
{
  // The first step is to load the data from the REMAP file
  ekat::ParameterList remap_params;
  remap_params.set("FILENAME",remap_filename);
  // We have 7 fields to load from the remap file:
  auto& remap_fields = remap_params.sublist("FIELDS");
  std::vector<std::string> remap_vars = {"S", "col", "row", "yc_a", "xc_a", "yc_b", "xc_b"};
  remap_fields.set("Number of Fields",remap_vars.size());
  int ifield = 0;
  for (auto f : remap_vars) {
    remap_fields.set(ekat::strint("field",ifield+1),f);
    ++ifield;
  }
  remap_params.print();
  // Now read remap data using these parameters.
}

} // namespace spa
} // namespace scream

#endif // SPA_HORIZONTAL_INTERPOLATION_IMPL_HPP
