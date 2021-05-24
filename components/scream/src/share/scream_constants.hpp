#ifndef SCREAM_CONSTANTS_HPP
#define SCREAM_CONSTANTS_HPP

#include <limits>

namespace scream {

/*
 * Mathematical and physical constants common across SCREAM.
 *
 */

template <typename Scalar>
struct Constants
{
  // Mathematical constants
  static constexpr Scalar zero          = 0.0;
  static constexpr Scalar one           = 1;
  static constexpr Scalar Pi            = 3.14159265358979;

  // Physical constants
  static constexpr Scalar gravity       = 9.80616;      // Gravity acceleration   [ m / s^2 ]
  static constexpr Scalar R             = 8.31446;      // Ideal gas constant [ J / (kg mol) ]

  // Atm-specific physical constants
  static constexpr Scalar P0            = 100000.0;   // Reference atm pressure [ Pa ]

  static constexpr Scalar M_dry_air     = 0.028966;   // Molar mass of dry air [ kg / mol ]
  static constexpr Scalar M_water_vapor = 18.016;     // Molar mass of water vapor [ g / mol ]

  static constexpr Scalar R_dry_air     = R / M_dry_air;      // Dry air gas constant [ J / (kg K) ]
  static constexpr Scalar R_water_vapor = R / M_water_vapor;  // Dry air gas constant [ J / (kg K) ]

  static constexpr Scalar cp_dry_air    = 1004.64;    // Dry air specific heat at constant pressure [ J/kg ]

  static constexpr Scalar rho_water     = 1000.0;     // Water density  [ kg ]
  static constexpr Scalar rho_ice       = 917.0;      // Ice density at 0 C from Wallace+Hobbes 1977 [ kg ]

  // Miscellanea constants
  static constexpr Scalar eps           = std::numeric_limits<Scalar>::epsilon(); // Machine precision
};

} // namespace scream

#endif // SCREAM_CONSTANTS_HPP
