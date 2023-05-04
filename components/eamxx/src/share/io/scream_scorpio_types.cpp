#include "scream_scorpio_types.hpp"

#include <ekat/ekat_assert.hpp>

namespace scream {
namespace scorpio {

std::string e2str (const FileMode mode)
{
  auto mode_int = static_cast<typename std::underlying_type<FileMode>::type>(mode);
  std::string s;
  switch (mode) {
    case Unset:  s = "UNSET";  break;
    case Read:   s = "READ";   break;
    case Write:  s = "WRITE";  break;
    case Append: s = "APPEND"; break;
    default:
      EKAT_ERROR_MSG (
          "Error! Unsupported/unrecognized FileMode value.\n"
          " - value: " + std::to_string(mode_int) + "\n");
  }
  return s;
}

} // namespace scorpio
} // namespace scream
