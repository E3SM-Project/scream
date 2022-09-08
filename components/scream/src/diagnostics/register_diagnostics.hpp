#ifndef SCREAM_REGISTER_DIAGNOSTICS_HPP
#define SCREAM_REGISTER_DIAGNOSTICS_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
// Include all diagnostics
#include "diagnostics/potential_temperature.hpp"
#include "diagnostics/atm_density.hpp"
#include "diagnostics/exner.hpp"
#include "diagnostics/virtual_temperature.hpp"
#include "diagnostics/vertical_layer_interface.hpp"
#include "diagnostics/vertical_layer_thickness.hpp"
#include "diagnostics/vertical_layer_midpoint.hpp"
#include "diagnostics/dry_static_energy.hpp"
#include "diagnostics/sea_level_pressure.hpp"
#include "diagnostics/liquid_water_path.hpp"
#include "diagnostics/ice_water_path.hpp"
#include "diagnostics/vapor_water_path.hpp"
#include "diagnostics/rain_water_path.hpp"

namespace scream {

inline void register_diagnostics () {
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  diag_factory.register_product("PotentialTemperature",&create_atmosphere_diagnostic<PotentialTemperatureDiagnostic>);
  diag_factory.register_product("AtmosphereDensity",&create_atmosphere_diagnostic<AtmDensityDiagnostic>);
  diag_factory.register_product("Exner",&create_atmosphere_diagnostic<ExnerDiagnostic>);
  diag_factory.register_product("VirtualTemperature",&create_atmosphere_diagnostic<VirtualTemperatureDiagnostic>);
  diag_factory.register_product("VerticalLayerInterface",&create_atmosphere_diagnostic<VerticalLayerInterfaceDiagnostic>);
  diag_factory.register_product("VerticalLayerThickness",&create_atmosphere_diagnostic<VerticalLayerThicknessDiagnostic>);
  diag_factory.register_product("VerticalLayerMidpoint",&create_atmosphere_diagnostic<VerticalLayerMidpointDiagnostic>);
  diag_factory.register_product("DryStaticEnergy",&create_atmosphere_diagnostic<DryStaticEnergyDiagnostic>);
  diag_factory.register_product("SeaLevelPressure",&create_atmosphere_diagnostic<SeaLevelPressureDiagnostic>);
  diag_factory.register_product("LiqWaterPath",&create_atmosphere_diagnostic<LiqWaterPathDiagnostic>);
  diag_factory.register_product("IceWaterPath",&create_atmosphere_diagnostic<IceWaterPathDiagnostic>);
  diag_factory.register_product("VapWaterPath",&create_atmosphere_diagnostic<VapWaterPathDiagnostic>);
  diag_factory.register_product("RainWaterPath",&create_atmosphere_diagnostic<RainWaterPathDiagnostic>);
}

} // namespace scream
#endif // SCREAM_REGISTER_DIAGNOSTICS_HPP
