# Volatility Basis Set (VBS) approach for SOA

## Overview

A modified volatility basis set (VBS) approach is used for both SOA precursor gases and particulate SOA that includes gas‐phase functionalization/fragmentation and particle‐phase oligomerization similar to FragNVSOA configuration of Shrivastava et al. (2015).[@shrivastava_global_2015] It includes a detailed treatment of SOA precursor gas chemistry including multigenerational aging via fragmentation and functionalization reactions, particle‐phase oligomerization that generates low “effective volatility” SOA, and particle‐phase loss by photolysis. The sources of SOA include natural biogenic, anthropogenic and biomass burning organic gases that are oxidized to form lower volatility species and undergo dynamic gas-particle partitioning to form SOA. Biogenic SOA is formed by oxidation of isoprene (ISOP_VBS) and monoterpene (C10H16) volatile organic compounds (VOCs). Emissions of anthropogenic and biomass burning organic gases in the range of intermediate volatility organics (IVOCs, referred to as SOAG0) are estimated as total primary organic matter (POM) emitted from these sources multiplied by specified tunable factors. IVOCs undergo multigenerational aging with OH radicals forming SOA corresponding to anthropogenic and biomass burning sources. Photolysis of SOA is included as a chemical sink in the particle phase, in addition to dry and wet removal sinks. The photolysis rate constant of particle-phase SOA is assumed to be 0.04% of typical NO2 photolysis frequencies following Hodzic et al. (2016).[@hodzic_rethinking_2016] The emissions of VBS SOA related gas species and oxidants (prescribed) read from a file are documented in the [User Guide](../user-guide/index.md).