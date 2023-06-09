#ifndef SHOC_MAIN_IMPL_HPP
#define SHOC_MAIN_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include "ekat/kokkos/ekat_subview_utils.hpp"

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc shoc_main. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
Int Functions<S,D>::shoc_init(
  const Int&                  nbot_shoc,
  const Int&                  ntop_shoc,
  const view_1d<const Spack>& pref_mid)
{
  // This function calculates the maximum number of levels
  // in pbl from surface

  using ExeSpace = typename KT::ExeSpace;
  view_1d<Int> npbl_d("npbl",1);

  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

    const Scalar pblmaxp = SC::pblmaxp;

    Int npbl_val = 1;

    const int begin_pack_indx = ntop_shoc/Spack::n;
    const int end_pack_indx   = nbot_shoc/Spack::n+1;
    Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, begin_pack_indx, end_pack_indx),
                                                    [&] (const Int& k, Int& local_max) {
      auto range = ekat::range<IntSmallPack>(k*Spack::n);
      auto condition = (range >= ntop_shoc && range < nbot_shoc);
      if (condition.any()) {
        condition = condition && pref_mid(k) >= pblmaxp;
      }

      auto levels_from_surface = nbot_shoc - range;
      levels_from_surface.set(!condition, 1);

      if (local_max < ekat::max(levels_from_surface))
        local_max = ekat::max(levels_from_surface);

    }, Kokkos::Max<Int>(npbl_val));

    npbl_d(0) = npbl_val;
  });

  const auto host_view = Kokkos::create_mirror_view(npbl_d);
  Kokkos::deep_copy(host_view, npbl_d);

  return host_view(0);
}

template<typename S, typename D>
Int Functions<S,D>::shoc_main(
  const Int&               shcol,               // Number of SHOC columns in the array
  const Int&               nlev,                // Number of levels
  const Int&               nlevi,               // Number of levels on interface grid
  const Int&               npbl,                // Maximum number of levels in pbl from surface
  const Int&               nadv,                // Number of times to loop SHOC
  const Int&               num_qtracers,        // Number of tracers
  const Scalar&            dtime,               // SHOC timestep [s]
  WorkspaceMgr&            workspace_mgr,       // WorkspaceManager for local variables
  const SHOCInput&         shoc_input,          // Input
  const SHOCInputOutput&   shoc_input_output,   // Input/Output
  const SHOCOutput&        shoc_output,         // Output
  const SHOCHistoryOutput& shoc_history_output  // Output (diagnostic)
#ifdef SCREAM_SMALL_KERNELS
  , const SHOCTemporaries& shoc_temporaries     // Temporaries for small kernels
#endif
                              )
{
  // Start timer
  auto start = std::chrono::steady_clock::now();

#ifndef SCREAM_SMALL_KERNELS
  using ExeSpace = typename KT::ExeSpace;

  // SHOC main loop
  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    const Scalar dx{shoc_input.dx(i)};
    const Scalar dy{shoc_input.dy(i)};
    const Scalar wthl_sfc{shoc_input.wthl_sfc(i)};
    const Scalar wqw_sfc{shoc_input.wqw_sfc(i)};
    const Scalar uw_sfc{shoc_input.uw_sfc(i)};
    const Scalar vw_sfc{shoc_input.vw_sfc(i)};
    const Scalar phis{shoc_input.phis(i)};
    Scalar& pblh = shoc_output.pblh(i);

    const auto zt_grid      = ekat::subview(shoc_input.zt_grid, i);
    const auto zi_grid      = ekat::subview(shoc_input.zi_grid, i);
    const auto pres         = ekat::subview(shoc_input.pres, i);
    const auto presi        = ekat::subview(shoc_input.presi, i);
    const auto pdel         = ekat::subview(shoc_input.pdel, i);
    const auto thv          = ekat::subview(shoc_input.thv, i);
    const auto w_field      = ekat::subview(shoc_input.w_field, i);
    const auto wtracer_sfc  = ekat::subview(shoc_input.wtracer_sfc, i);
    const auto inv_exner    = ekat::subview(shoc_input.inv_exner, i);
    const auto host_dse     = ekat::subview(shoc_input_output.host_dse, i);
    const auto tke          = ekat::subview(shoc_input_output.tke, i);
    const auto thetal       = ekat::subview(shoc_input_output.thetal, i);
    const auto qw           = ekat::subview(shoc_input_output.qw, i);
    const auto wthv_sec     = ekat::subview(shoc_input_output.wthv_sec, i);
    const auto tk           = ekat::subview(shoc_input_output.tk, i);
    const auto shoc_cldfrac = ekat::subview(shoc_input_output.shoc_cldfrac, i);
    const auto shoc_ql      = ekat::subview(shoc_input_output.shoc_ql, i);
    const auto shoc_ql2     = ekat::subview(shoc_output.shoc_ql2, i);
    const auto shoc_mix     = ekat::subview(shoc_history_output.shoc_mix, i);
    const auto w_sec        = ekat::subview(shoc_history_output.w_sec, i);
    const auto thl_sec      = ekat::subview(shoc_history_output.thl_sec, i);
    const auto qw_sec       = ekat::subview(shoc_history_output.qw_sec, i);
    const auto qwthl_sec    = ekat::subview(shoc_history_output.qwthl_sec, i);
    const auto wthl_sec     = ekat::subview(shoc_history_output.wthl_sec, i);
    const auto wqw_sec      = ekat::subview(shoc_history_output.wqw_sec, i);
    const auto wtke_sec     = ekat::subview(shoc_history_output.wtke_sec, i);
    const auto uw_sec       = ekat::subview(shoc_history_output.uw_sec, i);
    const auto vw_sec       = ekat::subview(shoc_history_output.vw_sec, i);
    const auto w3           = ekat::subview(shoc_history_output.w3, i);
    const auto wqls_sec     = ekat::subview(shoc_history_output.wqls_sec, i);
    const auto brunt        = ekat::subview(shoc_history_output.brunt, i);
    const auto isotropy     = ekat::subview(shoc_history_output.isotropy, i);

    const auto u_wind   = Kokkos::subview(shoc_input_output.horiz_wind, i, 0, Kokkos::ALL());
    const auto v_wind   = Kokkos::subview(shoc_input_output.horiz_wind, i, 1, Kokkos::ALL());
    const auto qtracers = Kokkos::subview(shoc_input_output.qtracers, i, Kokkos::ALL(), Kokkos::ALL());

    // Define temporary variables
    uview_1d<Spack> rho_zt, shoc_qv, dz_zt, dz_zi, tkh;
    workspace.template take_many_and_reset<5>(
      {"rho_zt", "shoc_qv", "dz_zt", "dz_zi", "tkh"},
      {&rho_zt, &shoc_qv, &dz_zt, &dz_zi, &tkh});

    // Local scalars
    Scalar se_b{0},   ke_b{0}, wv_b{0},   wl_b{0},
          se_a{0},   ke_a{0}, wv_a{0},   wl_a{0},
          ustar{0},  kbfs{0}, obklen{0}, ustar2{0}, wstar{0};

    // Scalarize some views for single entry access
    const auto s_thetal  = ekat::scalarize(thetal);
    const auto s_shoc_ql = ekat::scalarize(shoc_ql);
    const auto s_shoc_qv = ekat::scalarize(shoc_qv);

    // Compute integrals of static energy, kinetic energy, water vapor, and liquid water
    // for the computation of total energy before SHOC is called.  This is for an
    // effort to conserve energy since liquid water potential temperature (which SHOC
    // conserves) and static energy (which E3SM conserves) are not exactly equal.
    shoc_energy_integrals(team,nlev,host_dse,pdel,qw,shoc_ql,u_wind,v_wind, // Input
                          se_b,ke_b,wv_b,wl_b);                             // Output

    for (Int t=0; t<nadv; ++t) {
      // Check TKE to make sure values lie within acceptable
      // bounds after host model performs horizontal advection
      check_tke(team,nlev, // Input
                tke);      // Input/Output

      // Define vertical grid arrays needed for
      // vertical derivatives in SHOC, also
      // define air density (rho_zt)
      shoc_grid(team,nlev,nlevi,      // Input
                zt_grid,zi_grid,pdel, // Input
                dz_zt,dz_zi,rho_zt);  // Output

      // Compute the planetary boundary layer height, which is an
      // input needed for the length scale calculation.

      // Update SHOC water vapor,
      // to be used by the next two routines
      compute_shoc_vapor(team,nlev,qw,shoc_ql, // Input
                        shoc_qv);             // Output

      team.team_barrier();
      shoc_diag_obklen(uw_sfc,vw_sfc,     // Input
                      wthl_sfc, wqw_sfc, // Input
                      s_thetal(nlev-1),  // Input
                      s_shoc_ql(nlev-1), // Input
                      s_shoc_qv(nlev-1), // Input
                      ustar,kbfs,obklen); // Output

      pblintd(team,nlev,nlevi,npbl,     // Input
              zt_grid,zi_grid,thetal,   // Input
              shoc_ql,shoc_qv,u_wind,   // Input
              v_wind,ustar,obklen,kbfs, // Input
              shoc_cldfrac,             // Input
              workspace,                // Workspace
              pblh);                    // Output

      // Update the turbulent length scale
      shoc_length(team,nlev,nlevi,dx,dy, // Input
                  zt_grid,zi_grid,dz_zt, // Input
                  tke,thv,               // Input
                  workspace,             // Workspace
                  brunt,shoc_mix);       // Output

      // Advance the SGS TKE equation
      shoc_tke(team,nlev,nlevi,dtime,wthv_sec,    // Input
              shoc_mix,dz_zi,dz_zt,pres,u_wind,  // Input
              v_wind,brunt,obklen,zt_grid,       // Input
              zi_grid,pblh,                      // Input
              workspace,                         // Workspace
              tke,tk,tkh,                        // Input/Output
              isotropy);                         // Output

      // Update SHOC prognostic variables here
      // via implicit diffusion solver
      team.team_barrier();
      update_prognostics_implicit(team,nlev,nlevi,num_qtracers,dtime,dz_zt,   // Input
                                  dz_zi,rho_zt,zt_grid,zi_grid,tk,tkh,uw_sfc, // Input
                                  vw_sfc,wthl_sfc,wqw_sfc,wtracer_sfc,        // Input
                                  workspace,                                  // Workspace
                                  thetal,qw,qtracers,tke,u_wind,v_wind);   // Input/Output

      // Diagnose the second order moments
      diag_second_shoc_moments(team,nlev,nlevi,thetal,qw,u_wind,v_wind,   // Input
                              tke,isotropy,tkh,tk,dz_zi,zt_grid,zi_grid, // Input
                              shoc_mix,wthl_sfc,wqw_sfc,uw_sfc,vw_sfc,   // Input
                              ustar2,wstar,                              // Input/Output
                              workspace,                                 // Workspace
                              thl_sec,qw_sec,wthl_sec,wqw_sec,qwthl_sec, // Output
                              uw_sec,vw_sec,wtke_sec,w_sec);             // Output

      // Diagnose the third moment of vertical velocity,
      //  needed for the PDF closure
      diag_third_shoc_moments(team,nlev,nlevi,w_sec,thl_sec,wthl_sec, // Input
                              isotropy,brunt,thetal,tke,dz_zt,dz_zi,  // Input
                              zt_grid,zi_grid,                        // Input
                              workspace,                              // Workspace
                              w3);                                    // Output

      // Call the PDF to close on SGS cloud and turbulence
      team.team_barrier();
      shoc_assumed_pdf(team,nlev,nlevi,thetal,qw,w_field,thl_sec,qw_sec, // Input
                      wthl_sec,w_sec,wqw_sec,qwthl_sec,w3,pres,         // Input
                      zt_grid, zi_grid,                                 // Input
                      workspace,                                        // Workspace
                      shoc_cldfrac,shoc_ql,wqls_sec,wthv_sec,shoc_ql2); // Ouptut

      // Check TKE to make sure values lie within acceptable
      // bounds after vertical advection, etc.
      check_tke(team,nlev,tke);
    }

    // End SHOC parameterization

    // Use SHOC outputs to update the host model
    // temperature
    update_host_dse(team,nlev,thetal,shoc_ql, // Input
                    inv_exner,zt_grid,phis,   // Input
                    host_dse);                // Output

    team.team_barrier();
    shoc_energy_integrals(team,nlev,host_dse,pdel,  // Input
                          qw,shoc_ql,u_wind,v_wind, // Input
                          se_a,ke_a,wv_a,wl_a);     // Output

    shoc_energy_fixer(team,nlev,nlevi,dtime,nadv,zt_grid,zi_grid, // Input
                      se_b,ke_b,wv_b,wl_b,se_a,ke_a,wv_a,wl_a,    // Input
                      wthl_sfc,wqw_sfc,rho_zt,tke,presi,          // Input
                      workspace,                                  // Workspace
                      host_dse);                                  // Output

    // Remaining code is to diagnose certain quantities
    // related to PBL.  No answer changing subroutines
    // should be placed at this point onward.

    // Update PBLH, as other routines outside of SHOC
    // may require this variable.

    // Update SHOC water vapor, to be used by the next two routines
    compute_shoc_vapor(team,nlev,qw,shoc_ql, // Input
                      shoc_qv);             // Output

    team.team_barrier();
    shoc_diag_obklen(uw_sfc,vw_sfc,      // Input
                    wthl_sfc,wqw_sfc,   // Input
                    s_thetal(nlev-1),   // Input
                    s_shoc_ql(nlev-1),  // Input
                    s_shoc_qv(nlev-1),  // Input
                    ustar,kbfs,obklen); // Output

    pblintd(team,nlev,nlevi,npbl,zt_grid,   // Input
            zi_grid,thetal,shoc_ql,shoc_qv, // Input
            u_wind,v_wind,ustar,obklen,     // Input
            kbfs,shoc_cldfrac,              // Input
            workspace,                      // Workspace
            pblh);                          // Output

    // Release temporary variables from the workspace
    workspace.template release_many_contiguous<5>(
      {&rho_zt, &shoc_qv, &dz_zt, &dz_zi, &tkh});
  });
  Kokkos::fence();
#else
  const auto dx = shoc_input.dx;
  const auto dy = shoc_input.dy;
  const auto zt_grid = shoc_input.zt_grid;
  const auto zi_grid = shoc_input.zi_grid;
  const auto pres = shoc_input.pres;
  const auto presi =  shoc_input.presi;
  const auto pdel = shoc_input.pdel;
  const auto thv = shoc_input.thv;
  const auto w_field = shoc_input.w_field;
  const auto wthl_sfc = shoc_input.wthl_sfc;
  const auto wqw_sfc = shoc_input.wqw_sfc;
  const auto uw_sfc = shoc_input.uw_sfc;
  const auto vw_sfc = shoc_input.vw_sfc;
  const auto wtracer_sfc = shoc_input.wtracer_sfc;
  const auto inv_exner = shoc_input.inv_exner;
  const auto phis = shoc_input.phis;
  const auto host_dse = shoc_input_output.host_dse;
  const auto tke = shoc_input_output.tke;
  const auto thetal = shoc_input_output.thetal;
  const auto qw = shoc_input_output.qw;
  const auto wthv_sec = shoc_input_output.wthv_sec;
  const auto qtracers = shoc_input_output.qtracers;
  const auto tk = shoc_input_output.tk;
  const auto shoc_cldfrac = shoc_input_output.shoc_cldfrac;
  const auto pblh = shoc_output.pblh;
  const auto shoc_ql = shoc_input_output.shoc_ql;
  const auto shoc_ql2 = shoc_output.shoc_ql2;
  const auto shoc_mix = shoc_history_output.shoc_mix;
  const auto w_sec = shoc_history_output.w_sec;
  const auto thl_sec = shoc_history_output.thl_sec;
  const auto qw_sec = shoc_history_output.qw_sec;
  const auto qwthl_sec = shoc_history_output.qwthl_sec;
  const auto wthl_sec = shoc_history_output.wthl_sec;
  const auto wqw_sec = shoc_history_output.wqw_sec;
  const auto wtke_sec = shoc_history_output.wtke_sec;
  const auto uw_sec = shoc_history_output.uw_sec;
  const auto vw_sec = shoc_history_output.vw_sec;
  const auto w3 = shoc_history_output.w3;
  const auto wqls_sec = shoc_history_output.wqls_sec;
  const auto brunt = shoc_history_output.brunt;
  const auto isotropy = shoc_history_output.isotropy;
  const auto se_b = shoc_temporaries.se_b;
  const auto ke_b = shoc_temporaries.ke_b;
  const auto wv_b = shoc_temporaries.wv_b;
  const auto wl_b = shoc_temporaries.wl_b;
  const auto se_a = shoc_temporaries.se_a;
  const auto ke_a = shoc_temporaries.ke_a;
  const auto wv_a = shoc_temporaries.wv_a;
  const auto wl_a = shoc_temporaries.wl_a;
  const auto ustar = shoc_temporaries.ustar;
  const auto kbfs = shoc_temporaries.kbfs;
  const auto obklen = shoc_temporaries.obklen;
  const auto ustar2 = shoc_temporaries.ustar2;
  const auto wstar = shoc_temporaries.wstar;
  const auto rho_zt = shoc_temporaries.rho_zt;
  const auto shoc_qv = shoc_temporaries.shoc_qv;
  const auto dz_zt = shoc_temporaries.dz_zt;
  const auto dz_zi = shoc_temporaries.dz_zi;
  const auto tkh = shoc_temporaries.tkh;

  // Subview horizontal winds
  const auto u_wind = Kokkos::subview(shoc_input_output.horiz_wind, Kokkos::ALL(), 0, Kokkos::ALL());
  const auto v_wind = Kokkos::subview(shoc_input_output.horiz_wind, Kokkos::ALL(), 1, Kokkos::ALL());

  // Scalarize some views for single entry access
  const auto s_thetal  = ekat::scalarize(thetal);
  const auto s_shoc_ql = ekat::scalarize(shoc_ql);
  const auto s_shoc_qv = ekat::scalarize(shoc_qv);

  // Compute integrals of static energy, kinetic energy, water vapor, and liquid water
  // for the computation of total energy before SHOC is called.  This is for an
  // effort to conserve energy since liquid water potential temperature (which SHOC
  // conserves) and static energy (which E3SM conserves) are not exactly equal.
  shoc_energy_integrals_disp(shcol,nlev,host_dse,pdel,qw, // Input
                             shoc_ql,u_wind,v_wind,       // Input
                             se_b, ke_b, wv_b, wl_b);     // Output

  for (Int t=0; t<nadv; ++t) {
    // Check TKE to make sure values lie within acceptable
    // bounds after host model performs horizontal advection
    check_tke_disp(shcol,nlev, // Input
                   tke);       // Input/Output

    // Define vertical grid arrays needed for
    // vertical derivatives in SHOC, also
    // define air density (rho_zt)
    shoc_grid_disp(shcol,nlev,nlevi,     // Input
                   zt_grid,zi_grid,pdel, // Input
                   dz_zt,dz_zi,rho_zt);  // Output

    // Compute the planetary boundary layer height, which is an
    // input needed for the length scale calculation.

    // Update SHOC water vapor,
    // to be used by the next two routines
    compute_shoc_vapor_disp(shcol,nlev,qw,shoc_ql, // Input
                            shoc_qv);              // Output

    shoc_diag_obklen_disp(shcol, nlev,
                          uw_sfc,vw_sfc,      // Input
                          wthl_sfc, wqw_sfc,  // Input
                          s_thetal,           // Input
                          s_shoc_ql,          // Input
                          s_shoc_qv,          // Input
                          ustar,kbfs,obklen); // Output

    pblintd_disp(shcol,nlev,nlevi,npbl,    // Input
                 zt_grid,zi_grid,thetal,   // Input
                 shoc_ql,shoc_qv,u_wind,   // Input
                 v_wind,ustar,obklen,kbfs, // Input
                 shoc_cldfrac,             // Input
                 workspace_mgr,            // Workspace mgr
                 pblh);                    // Output

    // Update the turbulent length scale
    shoc_length_disp(shcol,nlev,nlevi,dx,dy, // Input
                     zt_grid,zi_grid,dz_zt,  // Input
                     tke,thv,                // Input
                     workspace_mgr,          // Workspace mgr
                     brunt,shoc_mix);        // Output

    // Advance the SGS TKE equation
    shoc_tke_disp(shcol,nlev,nlevi,dtime,wthv_sec,   // Input
                  shoc_mix,dz_zi,dz_zt,pres,u_wind,  // Input
                  v_wind,brunt,obklen,zt_grid,       // Input
                  zi_grid,pblh,                      // Input
                  workspace_mgr,                     // Workspace mgr
                  tke,tk,tkh,                        // Input/Output
                  isotropy);                         // Output

    // Update SHOC prognostic variables here
    // via implicit diffusion solver
    update_prognostics_implicit_disp(shcol,nlev,nlevi,num_qtracers,dtime,dz_zt,  // Input
                                     dz_zi,rho_zt,zt_grid,zi_grid,tk,tkh,uw_sfc, // Input
                                     vw_sfc,wthl_sfc,wqw_sfc,wtracer_sfc,        // Input
                                     workspace_mgr,                              // Workspace mgr
                                     thetal,qw,qtracers,tke,u_wind,v_wind);      // Input/Output

    // Diagnose the second order moments
    diag_second_shoc_moments_disp(shcol,nlev,nlevi,thetal,qw,u_wind,v_wind,  // Input
                                  tke,isotropy,tkh,tk,dz_zi,zt_grid,zi_grid, // Input
                                  shoc_mix,wthl_sfc,wqw_sfc,uw_sfc,vw_sfc,   // Input
                                  ustar2,wstar,                              // Input/Output
                                  workspace_mgr,                             // Workspace
                                  thl_sec,qw_sec,wthl_sec,wqw_sec,qwthl_sec, // Output
                                  uw_sec,vw_sec,wtke_sec,w_sec);             // Output

    // Diagnose the third moment of vertical velocity,
    //  needed for the PDF closure
    diag_third_shoc_moments_disp(shcol,nlev,nlevi,w_sec,thl_sec,wthl_sec, // Input
                                 isotropy,brunt,thetal,tke,dz_zt,dz_zi,   // Input
                                 zt_grid,zi_grid,                         // Input
                                 workspace_mgr,                           // Workspace mgr
                                 w3);                                     // Output

    // Call the PDF to close on SGS cloud and turbulence
    shoc_assumed_pdf_disp(shcol,nlev,nlevi,thetal,qw,w_field,thl_sec,qw_sec, // Input
                          wthl_sec,w_sec,wqw_sec,qwthl_sec,w3,pres,          // Input
                          zt_grid, zi_grid,                                  // Input
                          workspace_mgr,                                     // Workspace mgr
                          shoc_cldfrac,shoc_ql,wqls_sec,wthv_sec,shoc_ql2);  // Ouptut

    // Check TKE to make sure values lie within acceptable
    // bounds after vertical advection, etc.
    check_tke_disp(shcol,nlev,tke);
  }

  // End SHOC parameterization

  // Use SHOC outputs to update the host model
  // temperature
  update_host_dse_disp(shcol,nlev,thetal,shoc_ql, // Input
                       inv_exner,zt_grid,phis,    // Input
                       host_dse);                 // Output

  shoc_energy_integrals_disp(shcol,nlev,host_dse,pdel,  // Input
                             qw,shoc_ql,u_wind,v_wind,  // Input
                             se_a,ke_a,wv_a,wl_a);      // Output

  shoc_energy_fixer_disp(shcol,nlev,nlevi,dtime,nadv,zt_grid,zi_grid, // Input
                         se_b,ke_b,wv_b,wl_b,se_a,ke_a,wv_a,wl_a,     // Input
                         wthl_sfc,wqw_sfc,rho_zt,tke,presi,           // Input
                         workspace_mgr,                               // Workspace
                         host_dse);                                   // Output

  // Remaining code is to diagnose certain quantities
  // related to PBL.  No answer changing subroutines
  // should be placed at this point onward.

  // Update PBLH, as other routines outside of SHOC
  // may require this variable.

  // Update SHOC water vapor, to be used by the next two routines
  compute_shoc_vapor_disp(shcol,nlev,qw,shoc_ql, // Input
                          shoc_qv);              // Output

  shoc_diag_obklen_disp(shcol, nlev, uw_sfc,vw_sfc, // Input
                        wthl_sfc,wqw_sfc,           // Input
                        s_thetal,                   // Input
                        s_shoc_ql,                  // Input
                        s_shoc_qv,                  // Input
                        ustar,kbfs,obklen);         // Output

  pblintd_disp(shcol,nlev,nlevi,npbl,zt_grid,  // Input
               zi_grid,thetal,shoc_ql,shoc_qv, // Input
               u_wind,v_wind,ustar,obklen,     // Input
               kbfs,shoc_cldfrac,              // Input
               workspace_mgr,                  // Workspace mgr
               pblh);                          // Output
#endif

  auto finish = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
  return duration.count();
}

} // namespace shoc
} // namespace scream

#endif
