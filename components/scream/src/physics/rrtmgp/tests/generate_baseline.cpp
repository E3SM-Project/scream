#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/rrtmgp/mo_garand_atmos_io.h"
#include "physics/rrtmgp/rrtmgp_test_utils.hpp"
#include "share/scream_types.hpp"
#include "share/scream_session.hpp"

#include "cpp/rrtmgp/mo_gas_concentrations.h"

#include "YAKL.h"

#include <iostream>
#include <cmath>

using namespace scream;

int main (int argc, char** argv) {

    // Get filenames from command line
    if (argc != 3) {
        std::cout <<
            argv[0] << " [options] inputfile baseline\n"
            "Options:\n"
            "  (there are no options)\n";
        return 1;
    }
    std::string inputfile(argv[argc-2]);
    std::string baseline(argv[argc-1]);

    // Initialize yakl
    yakl::init();

    // Get reference fluxes from input file; do this here so we can get ncol dimension
    real2d sw_flux_up_ref;
    real2d sw_flux_dn_ref;
    real2d sw_flux_dn_dir_ref;
    real2d lw_flux_up_ref;
    real2d lw_flux_dn_ref;
    std::cout << "read_fluxes..." << std::endl;
    rrtmgpTest::read_fluxes(inputfile, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dn_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

    // Get dimension sizes
    int ncol = sw_flux_up_ref.dimension[0];
    int nlev = sw_flux_up_ref.dimension[1];
    int nlay = nlev - 1;

    // Read in dummy Garand atmosphere; if this were an actual model simulation,
    // these would be passed as inputs to the driver
    // NOTE: set ncol to size of col_flx dimension in the input file. This is so
    // that we can compare to the reference data provided in that file. Note that
    // this will copy the first column of the input data (the first profile) ncol
    // times. We will then fill some fraction of these columns with clouds for
    // the test problem.
    real2d p_lay ("p_lay", ncol, nlay);
    real2d t_lay ("t_lay", ncol, nlay);
    real2d p_lev ("p_lev", ncol, nlay+1);
    real2d t_lev ("t_lev", ncol, nlay+1);
    GasConcs gas_concs;
    read_atmos(inputfile, p_lay, t_lay, p_lev, t_lev, gas_concs, ncol);

    // Initialize the RRTMGP interface; this will read in the k-distribution
    // data that contains information about absorption coefficients for gases
    std::cout << "rrtmgp_initialize..." << std::endl;
    rrtmgp::rrtmgp_initialize(gas_concs);

    // Setup dummy all-sky problem
    real1d sfc_alb_dir_vis ("sfc_alb_dir_vis", ncol);
    real1d sfc_alb_dir_nir ("sfc_alb_dir_nir", ncol);
    real1d sfc_alb_dif_vis ("sfc_alb_dif_vis", ncol);
    real1d sfc_alb_dif_nir ("sfc_alb_dif_nir", ncol);
    real1d mu0 ("mu0", ncol);
    real2d lwp ("lwp", ncol, nlay);
    real2d iwp ("iwp", ncol, nlay);
    real2d rel ("rel", ncol, nlay);
    real2d rei ("rei", ncol, nlay);
    real2d cld ("cld", ncol, nlay);
    rrtmgpTest::dummy_atmos(
        inputfile, ncol, p_lay, t_lay,
        sfc_alb_dir_vis, sfc_alb_dir_nir,
        sfc_alb_dif_vis, sfc_alb_dif_nir,
        mu0,
        lwp, iwp, rel, rei, cld
    );

    // Setup flux outputs; In a real model run, the fluxes would be
    // input/outputs into the driver (persisting between calls), and
    // we would just have to setup the pointers to them in the
    // FluxesBroadband object
    real2d sw_flux_up ("sw_flux_up" ,ncol,nlay+1);
    real2d sw_flux_dn ("sw_flux_dn" ,ncol,nlay+1);
    real2d sw_flux_dn_dir("sw_flux_dn_dir",ncol,nlay+1);
    real2d lw_flux_up ("lw_flux_up" ,ncol,nlay+1);
    real2d lw_flux_dn ("lw_flux_dn" ,ncol,nlay+1);

    // Compute band-by-band surface_albedos.
    const auto nswbands = 14;
    real2d sfc_alb_dir("sfc_alb_dir", ncol, nswbands);
    real2d sfc_alb_dif("sfc_alb_dif", ncol, nswbands);
    rrtmgp::compute_band_by_band_surface_albedos(
      ncol, nswbands,
      sfc_alb_dir_vis, sfc_alb_dir_nir,
      sfc_alb_dif_vis, sfc_alb_dif_nir,
      sfc_alb_dir, sfc_alb_dif);

    // Run RRTMGP standalone codes and compare with AD run
    // Do something interesting here...
    // NOTE: these will get replaced with AD stuff that handles these
    std::cout << "rrtmgp_main..." << std::endl;
    rrtmgp::rrtmgp_main(
        ncol, nlay,
        p_lay, t_lay, p_lev, t_lev, gas_concs,
        sfc_alb_dir, sfc_alb_dif, mu0,
        lwp, iwp, rel, rei,
        sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
        lw_flux_up, lw_flux_dn
    );

    // Write fluxes
    std::cout << "write_fluxes..." << std::endl;
    rrtmgpTest::write_fluxes(
        baseline, 
        sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
        lw_flux_up, lw_flux_dn
    );

    // Clean up from test; this is probably not necessary, these things
    // should be deallocated when they fall out of scope, but we should be
    // good citizens and clean up our mess.
    p_lay.deallocate();
    t_lay.deallocate();
    p_lev.deallocate();
    t_lev.deallocate();
    sfc_alb_dir_vis.deallocate();
    sfc_alb_dir_nir.deallocate();
    sfc_alb_dif_vis.deallocate();
    sfc_alb_dif_nir.deallocate();
    sfc_alb_dir.deallocate();
    sfc_alb_dif.deallocate();
    mu0.deallocate();
    lwp.deallocate();
    iwp.deallocate();
    rel.deallocate();
    rei.deallocate();
    cld.deallocate();
    sw_flux_up_ref.deallocate();
    sw_flux_dn_ref.deallocate();
    sw_flux_dn_dir_ref.deallocate();
    lw_flux_up_ref.deallocate();
    lw_flux_dn_ref.deallocate();
    sw_flux_up.deallocate();
    sw_flux_dn.deallocate();
    sw_flux_dn_dir.deallocate();
    lw_flux_up.deallocate();
    lw_flux_dn.deallocate();

    gas_concs.reset();
    rrtmgp::rrtmgp_finalize();
    yakl::finalize();

    return 0;
}
