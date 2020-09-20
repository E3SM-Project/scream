module shoc_iso_c
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from scream c++ to shoc fortran.
!

contains
  subroutine append_precision(string, prefix)

    character(kind=c_char, len=128), intent(inout) :: string
    character(*), intent(in) :: prefix
    real(kind=c_real) :: s

    write (string, '(a,i1,a1)') prefix, sizeof(s), C_NULL_CHAR
  end subroutine append_precision

  subroutine shoc_init_c(nlev, gravit, rair, rh2o, cpair, &
                         zvir, latvap, latice, karman) bind(c)
    use shoc, only: shoc_init, npbl

    integer(kind=c_int), value, intent(in) :: nlev ! number of levels

    real(kind=c_real), value, intent(in)  :: gravit ! gravity
    real(kind=c_real), value, intent(in)  :: rair   ! dry air gas constant
    real(kind=c_real), value, intent(in)  :: rh2o   ! water vapor gas constant
    real(kind=c_real), value, intent(in)  :: cpair  ! specific heat of dry air
    real(kind=c_real), value, intent(in)  :: zvir   ! rh2o/rair - 1
    real(kind=c_real), value, intent(in)  :: latvap ! latent heat of vaporization
    real(kind=c_real), value, intent(in)  :: latice ! latent heat of fusion
    real(kind=c_real), value, intent(in)  :: karman ! Von Karman's constant

    real(kind=c_real) :: pref_mid(nlev) ! unused values

    pref_mid = 0
    call shoc_init(nlev, gravit, rair, rh2o, cpair, &
                   zvir, latvap, latice, karman, &
                   pref_mid, nlev, 1)
    npbl = nlev ! set pbl layer explicitly so we don't need pref_mid.
  end subroutine shoc_init_c

  subroutine shoc_main_c(shcol,nlev,nlevi,dtime,nadv,host_dx, host_dy, thv,  &
     zt_grid, zi_grid, pres, presi, pdel, wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &
     wtracer_sfc, num_qtracers, w_field, exner,phis, host_dse, tke, thetal,  &
     qw, u_wind, v_wind, qtracers, wthv_sec, tkh, tk, shoc_ql, shoc_cldfrac, &
     pblh, shoc_mix, isotropy, w_sec, thl_sec, qw_sec, qwthl_sec, wthl_sec,  &
     wqw_sec, wtke_sec, uw_sec, vw_sec, w3, wqls_sec, brunt, shoc_ql2) bind(C)
    use shoc, only : shoc_main

    integer(kind=c_int), value, intent(in) :: shcol, nlev, nlevi, num_qtracers, nadv
    real(kind=c_real), value, intent(in) :: dtime
    real(kind=c_real), intent(in), dimension(shcol) :: host_dx, host_dy
    real(kind=c_real), intent(in), dimension(shcol, nlev) :: zt_grid
    real(kind=c_real), intent(in), dimension(shcol, nlevi) :: zi_grid
    real(kind=c_real), intent(in), dimension(shcol, nlev) :: pres
    real(kind=c_real), intent(in), dimension(shcol, nlevi) :: presi
    real(kind=c_real), intent(in), dimension(shcol, nlev) :: pdel, thv, w_field
    real(kind=c_real), intent(in), dimension(shcol) :: wthl_sfc, wqw_sfc, uw_sfc, vw_sfc
    real(kind=c_real), intent(in), dimension(shcol, num_qtracers) :: wtracer_sfc
    real(kind=c_real), intent(in), dimension(shcol, nlev) :: exner
    real(kind=c_real), intent(in), dimension(shcol) :: phis

    real(kind=c_real), intent(inout), dimension(shcol, nlev) :: host_dse, tke, &
       thetal, qw, u_wind, v_wind, wthv_sec
    real(kind=c_real), intent(inout), dimension(shcol, nlev, num_qtracers) :: qtracers
    real(kind=c_real), intent(inout), dimension(shcol, nlev) :: tk, tkh

    real(kind=c_real), intent(out), dimension(shcol, nlev) :: shoc_cldfrac, shoc_ql
    real(kind=c_real), intent(out), dimension(shcol) :: pblh
    real(kind=c_real), intent(out), dimension(shcol, nlev) :: shoc_mix, w_sec
    real(kind=c_real), intent(out), dimension(shcol, nlevi) :: thl_sec, qw_sec, &
       qwthl_sec, wthl_sec, wqw_sec, wtke_sec, uw_sec, vw_sec, w3
    real(kind=c_real), intent(out), dimension(shcol, nlev) :: wqls_sec, isotropy, &
         brunt,shoc_ql2

    call shoc_main(shcol, nlev, nlevi, dtime, nadv, host_dx, host_dy, thv,   &
     zt_grid, zi_grid, pres, presi, pdel, wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &
     wtracer_sfc, num_qtracers, w_field, exner, phis, host_dse, tke, thetal, &
     qw, u_wind, v_wind, qtracers, wthv_sec, tkh, tk, shoc_ql, shoc_cldfrac, &
     pblh, shoc_mix, isotropy, w_sec, thl_sec, qw_sec, qwthl_sec, wthl_sec,  &
     wqw_sec, wtke_sec, uw_sec, vw_sec, w3, wqls_sec, brunt,shoc_ql2 )
  end subroutine shoc_main_c

  subroutine shoc_use_cxx_c(arg_use_cxx) bind(C)
    use shoc, only: use_cxx

    logical(kind=c_bool), value, intent(in) :: arg_use_cxx

    use_cxx = arg_use_cxx
  end subroutine shoc_use_cxx_c

  subroutine shoc_grid_c(shcol,nlev,nlevi,zt_grid,zi_grid,pdel,dz_zt,dz_zi,rho_zt) bind (C)
    use shoc, only: shoc_grid

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: nlevi
    real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
    real(kind=c_real), intent(in) :: zi_grid(shcol,nlevi)
    real(kind=c_real), intent(in) :: pdel(shcol,nlev)

    real(kind=c_real), intent(out) :: dz_zt(shcol,nlev)
    real(kind=c_real), intent(out) :: dz_zi(shcol,nlevi)
    real(kind=c_real), intent(out) :: rho_zt(shcol,nlev)

    call shoc_grid(shcol,nlev,nlevi,zt_grid,zi_grid,pdel,dz_zt,dz_zi,rho_zt)

  end subroutine shoc_grid_c

  subroutine calc_shoc_varorcovar_c(&
       shcol,nlev,nlevi,tunefac,&                ! Input
       isotropy_zi,tkh_zi,dz_zi,invar1,invar2,&  ! Input
       varorcovar)bind (C)                       ! Input/Output

      use shoc, only: calc_shoc_varorcovar

      integer(kind=c_int), intent(in), value :: shcol
      integer(kind=c_int), intent(in), value :: nlev
      integer(kind=c_int), intent(in), value :: nlevi
      real(kind=c_real), intent(in), value :: tunefac
      real(kind=c_real), intent(in) :: isotropy_zi(shcol,nlevi)
      real(kind=c_real), intent(in) :: tkh_zi(shcol,nlevi)
      real(kind=c_real), intent(in) :: dz_zi(shcol,nlevi)
      real(kind=c_real), intent(in) :: invar1(shcol,nlev)
      real(kind=c_real), intent(in) :: invar2(shcol,nlev)

      real(kind=c_real), intent(inout) :: varorcovar(shcol,nlevi)

      call calc_shoc_varorcovar(&
           shcol,nlev,nlevi,tunefac,&               ! Input
           isotropy_zi,tkh_zi,dz_zi,invar1,invar2,& ! Input
           varorcovar)

  end subroutine calc_shoc_varorcovar_c

  subroutine integ_column_stability_c(nlev, shcol, dz_zt, pres, brunt, brunt_int) bind (C)
    use shoc, only: integ_column_stability

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: dz_zt(shcol,nlev)
    real(kind=c_real), intent(in) :: pres(shcol,nlev)
    real(kind=c_real), intent(in) :: brunt(shcol,nlev)

    real(kind=c_real), intent(out) :: brunt_int(shcol)

    call integ_column_stability(nlev, shcol, dz_zt, pres, brunt, brunt_int)

  end subroutine integ_column_stability_c

  subroutine compute_shr_prod_c(nlevi, nlev, shcol, dz_zi, u_wind, v_wind, sterm) bind (C)
    use shoc, only: compute_shr_prod

    integer(kind=c_int), intent(in), value :: nlevi
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: dz_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: u_wind(shcol,nlev)
    real(kind=c_real), intent(in) :: v_wind(shcol,nlev)

    real(kind=c_real), intent(out) :: sterm(shcol,nlevi)

    call compute_shr_prod(nlevi, nlev, shcol, dz_zi, u_wind, v_wind, sterm)

  end subroutine compute_shr_prod_c

  subroutine isotropic_ts_c(nlev, shcol, brunt_int, tke, a_diss, brunt, isotropy) bind (C)
    use shoc, only: isotropic_ts

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: brunt_int(shcol)
    real(kind=c_real), intent(in) :: tke(shcol,nlev)
    real(kind=c_real), intent(in) :: a_diss(shcol,nlev)
    real(kind=c_real), intent(in) :: brunt(shcol,nlev)

    real(kind=c_real), intent(out) :: isotropy(shcol,nlev)

    call isotropic_ts(nlev, shcol, brunt_int, tke, a_diss, brunt, isotropy)

  end subroutine isotropic_ts_c

  subroutine adv_sgs_tke_c(nlev, shcol, dtime, shoc_mix, wthv_sec, &
                           sterm_zt, tk, tke, a_diss) bind (C)
    use shoc, only: adv_sgs_tke

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in), value :: dtime
    real(kind=c_real), intent(in) :: shoc_mix(shcol,nlev)
    real(kind=c_real), intent(in) :: wthv_sec(shcol,nlev)
    real(kind=c_real), intent(in) :: sterm_zt(shcol,nlev)
    real(kind=c_real), intent(in) :: tk(shcol,nlev)

    real(kind=c_real), intent(inout) :: tke(shcol,nlev)

    real(kind=c_real), intent(out) :: a_diss(shcol,nlev)

    call adv_sgs_tke(nlev, shcol, dtime, shoc_mix, wthv_sec, &
                     sterm_zt, tk, tke, a_diss)

  end subroutine adv_sgs_tke_c

  subroutine eddy_diffusivities_c(nlev, shcol, obklen, pblh, zt_grid, &
                          shoc_mix, sterm_zt, isotropy, tke, tkh, tk) bind (C)
    use shoc, only: eddy_diffusivities

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: obklen(shcol)
    real(kind=c_real), intent(in) :: pblh(shcol)
    real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
    real(kind=c_real), intent(in) :: shoc_mix(shcol,nlev)
    real(kind=c_real), intent(in) :: sterm_zt(shcol,nlev)
    real(kind=c_real), intent(in) :: isotropy(shcol,nlev)
    real(kind=c_real), intent(in) :: tke(shcol,nlev)

    real(kind=c_real), intent(out) :: tkh(shcol,nlev)
    real(kind=c_real), intent(out) :: tk(shcol,nlev)

    call eddy_diffusivities(nlev, shcol, obklen, pblh, zt_grid, &
     shoc_mix, sterm_zt, isotropy, tke, tkh, tk)

  end subroutine eddy_diffusivities_c

  subroutine update_host_dse_c(shcol, nlev, thlm, shoc_ql, exner, zt_grid, &
                               phis, host_dse) bind (C)
    use shoc, only: update_host_dse

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    real(kind=c_real), intent(in) :: thlm(shcol,nlev)
    real(kind=c_real), intent(in) :: shoc_ql(shcol,nlev)
    real(kind=c_real), intent(in) :: exner(shcol,nlev)
    real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
    real(kind=c_real), intent(in) :: phis(shcol)

    real(kind=c_real), intent(out) :: host_dse(shcol,nlev)

    call update_host_dse(shcol, nlev, thlm, shoc_ql, exner, zt_grid, &
                           phis, host_dse)

  end subroutine update_host_dse_c

  subroutine shoc_energy_integrals_c(shcol, nlev, host_dse, pdel,&
                                     rtm, rcm, u_wind, v_wind,&
                                     se_int, ke_int, wv_int, wl_int) bind (C)
    use shoc, only: shoc_energy_integrals

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    real(kind=c_real), intent(in) :: host_dse(shcol,nlev)
    real(kind=c_real), intent(in) :: pdel(shcol,nlev)
    real(kind=c_real), intent(in) :: rtm(shcol,nlev)
    real(kind=c_real), intent(in) :: rcm(shcol,nlev)
    real(kind=c_real), intent(in) :: u_wind(shcol,nlev)
    real(kind=c_real), intent(in) :: v_wind(shcol,nlev)

    real(kind=c_real), intent(out) :: se_int(shcol)
    real(kind=c_real), intent(out) :: ke_int(shcol)
    real(kind=c_real), intent(out) :: wv_int(shcol)
    real(kind=c_real), intent(out) :: wl_int(shcol)

    call shoc_energy_integrals(shcol, nlev, host_dse, pdel,&
                               rtm, rcm, u_wind, v_wind,&
                               se_int, ke_int, wv_int, wl_int)

  end subroutine shoc_energy_integrals_c
  
  subroutine shoc_energy_total_fixer_c(&
                                shcol,nlev,nlevi,dtime,nadv,&
                                zt_grid,zi_grid,&
                                se_b,ke_b,wv_b,wl_b,&
                                se_a,ke_a,wv_a,wl_a,&
                                wthl_sfc,wqw_sfc,rho_zt,&
                                te_a,te_b) bind (C)
    use shoc, only: shoc_energy_total_fixer

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: nlevi
    integer(kind=c_int), intent(in), value :: nadv
    real(kind=c_real), intent(in), value :: dtime
    
    real(kind=c_real), intent(in) :: se_b(shcol)
    real(kind=c_real), intent(in) :: ke_b(shcol)
    real(kind=c_real), intent(in) :: wv_b(shcol)
    real(kind=c_real), intent(in) :: wl_b(shcol)
    real(kind=c_real), intent(in) :: se_a(shcol)
    real(kind=c_real), intent(in) :: ke_a(shcol)
    real(kind=c_real), intent(in) :: wv_a(shcol)
    real(kind=c_real), intent(in) :: wl_a(shcol)
    real(kind=c_real), intent(in) :: wthl_sfc(shcol)
    real(kind=c_real), intent(in) :: wqw_sfc(shcol)
    real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
    real(kind=c_real), intent(in) :: zi_grid(shcol,nlevi)
    real(kind=c_real), intent(in) :: rho_zt(shcol,nlev)

    real(kind=c_real), intent(out) :: te_a(shcol)
    real(kind=c_real), intent(out) :: te_b(shcol)

    call shoc_energy_total_fixer(shcol,nlev,nlevi,dtime,nadv,&
                                 zt_grid,zi_grid,&
                                 se_b,ke_b,wv_b,wl_b,&
                                 se_a,ke_a,wv_a,wl_a,&
                                 wthl_sfc,wqw_sfc,rho_zt,&
                                 te_a,te_b)

  end subroutine shoc_energy_total_fixer_c

  subroutine shoc_energy_threshold_fixer_c(shcol,nlev,nlevi, &
                                     pint,tke,te_a,te_b, &
                                     se_dis,shoctop) bind (C)
    use shoc, only: shoc_energy_threshold_fixer

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: nlevi
    
    real(kind=c_real), intent(in) :: pint(shcol,nlevi)
    real(kind=c_real), intent(in) :: tke(shcol,nlev)
    real(kind=c_real), intent(in) :: te_a(shcol)
    real(kind=c_real), intent(in) :: te_b(shcol)
    
    real(kind=c_real), intent(out) :: se_dis(shcol)
    integer(kind=c_int), intent(out) :: shoctop(shcol)

    call shoc_energy_threshold_fixer(shcol,nlev,nlevi, &
                                     pint,tke,te_a,te_b, &
                                     se_dis,shoctop)

  end subroutine shoc_energy_threshold_fixer_c
  
  subroutine shoc_energy_dse_fixer_c(shcol,nlev, &
                                     se_dis,shoctop, &
                                     host_dse) bind (C)
    use shoc, only: shoc_energy_dse_fixer

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in) :: shoctop(shcol)
    real(kind=c_real), intent(in) :: se_dis(shcol)

    real(kind=c_real), intent(inout) :: host_dse(shcol,nlev)

    call shoc_energy_dse_fixer(shcol,nlev, &
                               se_dis,shoctop, &
                               host_dse)

  end subroutine shoc_energy_dse_fixer_c  

  subroutine calc_shoc_vertflux_c(shcol, nlev, nlevi, tkh_zi, dz_zi, invar, vertflux) bind (C)
    use shoc, only: calc_shoc_vertflux

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: nlevi
    real(kind=c_real), intent(in) :: tkh_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: dz_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: invar(shcol,nlev)

    real(kind=c_real), intent(inout) :: vertflux(shcol,nlevi)

    call calc_shoc_vertflux(shcol, nlev, nlevi, tkh_zi, dz_zi, invar, vertflux)

  end subroutine calc_shoc_vertflux_c

  subroutine compute_brunt_shoc_length_c(nlev,nlevi,shcol,dz_zt,thv,thv_zi,brunt) bind (C)
    use shoc, only: compute_brunt_shoc_length

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: nlevi
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: dz_zt(shcol,nlev)
    real(kind=c_real), intent(in) :: thv(shcol,nlev)  
    real(kind=c_real), intent(in) :: thv_zi(shcol,nlevi)

    real(kind=c_real), intent(out) :: brunt(shcol,nlev)

    call compute_brunt_shoc_length(nlev,nlevi,shcol,dz_zt,thv,thv_zi,brunt)

  end subroutine compute_brunt_shoc_length_c  
  
  
  subroutine compute_l_inf_shoc_length_c(nlev,shcol,zt_grid,dz_zt,tke,l_inf) bind (C)
    use shoc, only: compute_l_inf_shoc_length

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
    real(kind=c_real), intent(in) :: dz_zt(shcol,nlev)
    real(kind=c_real), intent(in) :: tke(shcol,nlev)  

    real(kind=c_real), intent(out) :: l_inf(shcol)

    call compute_l_inf_shoc_length(nlev,shcol,zt_grid,dz_zt,tke,l_inf)

  end subroutine compute_l_inf_shoc_length_c  
  
  
  subroutine compute_conv_vel_shoc_length_c(nlev,shcol,pblh,zt_grid,dz_zt,&
                                            thv,wthv_sec,conv_vel) bind (C)
    use shoc, only: compute_conv_vel_shoc_length

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: pblh(shcol)
    real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
    real(kind=c_real), intent(in) :: dz_zt(shcol,nlev)
    real(kind=c_real), intent(in) :: thv(shcol,nlev)
    real(kind=c_real), intent(in) :: wthv_sec(shcol,nlev)  

    real(kind=c_real), intent(out) :: conv_vel(shcol)

    call compute_conv_vel_shoc_length(nlev,shcol,pblh,zt_grid,dz_zt,thv,wthv_sec,conv_vel)

  end subroutine compute_conv_vel_shoc_length_c 
  
  subroutine compute_conv_time_shoc_length_c(shcol,pblh,conv_vel,tscale) bind (C)
    use shoc, only: compute_conv_time_shoc_length

    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: pblh(shcol)  
    real(kind=c_real), intent(inout) :: conv_vel(shcol)
    real(kind=c_real), intent(out) :: tscale(shcol)

    call compute_conv_time_shoc_length(shcol,pblh,conv_vel,tscale)

  end subroutine compute_conv_time_shoc_length_c 
  
  subroutine compute_shoc_mix_shoc_length_c(nlev,shcol,tke,brunt,tscale,&
                                         zt_grid,l_inf,shoc_mix) bind (C)
    use shoc, only: compute_shoc_mix_shoc_length

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: tke(shcol,nlev)
    real(kind=c_real), intent(in) :: brunt(shcol,nlev)
    real(kind=c_real), intent(in) :: tscale(shcol)
    real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
    real(kind=c_real), intent(in) :: l_inf(shcol) 

    real(kind=c_real), intent(out) :: shoc_mix(shcol,nlev)

    call compute_shoc_mix_shoc_length(nlev,shcol,tke,brunt,tscale,zt_grid,&
                                      l_inf,shoc_mix)

  end subroutine compute_shoc_mix_shoc_length_c 
  
  subroutine check_length_scale_shoc_length_c(nlev,shcol,host_dx,host_dy,shoc_mix) bind (C)
    use shoc, only: check_length_scale_shoc_length

    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: shcol
    real(kind=c_real), intent(in) :: host_dx(shcol)
    real(kind=c_real), intent(in) :: host_dy(shcol)  

    real(kind=c_real), intent(inout) :: shoc_mix(shcol,nlev)

    call check_length_scale_shoc_length(nlev,shcol,host_dx,host_dy,shoc_mix)

  end subroutine check_length_scale_shoc_length_c         

  subroutine shoc_diag_second_moments_srf_c(shcol, wthl, uw, vw, ustar2, wstar) bind(C)
   use shoc, only: diag_second_moments_srf

   ! argmens
   integer(kind=c_int), value, intent(in) :: shcol
   real(kind=c_real), intent(in)  :: wthl(shcol), uw(shcol), vw(shcol)
   real(kind=c_real), intent(out) :: ustar2(shcol), wstar(shcol)

   call diag_second_moments_srf(shcol, wthl, uw, vw, ustar2, wstar)
 end subroutine shoc_diag_second_moments_srf_c
 
  subroutine linear_interp_c(x1,x2,y1,y2,km1,km2,ncol,minthresh) bind (C)
    use shoc, only: linear_interp

    integer(kind=c_int), intent(in), value :: ncol
    integer(kind=c_int), intent(in), value :: km1
    integer(kind=c_int), intent(in), value :: km2
    real(kind=c_real), intent(in), value :: minthresh
    real(kind=c_real), intent(in) :: x1(ncol,km1)
    real(kind=c_real), intent(in) :: y1(ncol,km1)
    real(kind=c_real), intent(in) :: x2(ncol,km2)
    
    real(kind=c_real), intent(out) :: y2(ncol,km2)

    call linear_interp(x1,x2,y1,y2,km1,km2,ncol,minthresh)

  end subroutine linear_interp_c  
  
  subroutine shoc_assumed_pdf_tilda_to_real_c(w_first, sqrtw2, w1) bind (C)
    use shoc, only: shoc_assumed_pdf_tilda_to_real

    real(kind=c_real), intent(in), value :: w_first
    real(kind=c_real), intent(in), value :: sqrtw2
    
    real(kind=c_real), intent(inout) :: w1

    call shoc_assumed_pdf_tilda_to_real(w_first, sqrtw2, w1)

  end subroutine shoc_assumed_pdf_tilda_to_real_c   
  
  subroutine shoc_assumed_pdf_vv_parameters_c(w_first,w_sec,w3var,&         
                                              Skew_w,w1_1,w1_2,w2_1,w2_2,a) bind (C)
    use shoc, only: shoc_assumed_pdf_vv_parameters

    real(kind=c_real), intent(in), value :: w_first
    real(kind=c_real), intent(in), value :: w_sec
    real(kind=c_real), intent(in), value :: w3var
    
    real(kind=c_real), intent(out) :: Skew_w
    real(kind=c_real), intent(out) :: w1_1
    real(kind=c_real), intent(out) :: w1_2
    real(kind=c_real), intent(out) :: w2_1
    real(kind=c_real), intent(out) :: w2_2
    real(kind=c_real), intent(out) :: a

    call shoc_assumed_pdf_vv_parameters(w_first,w_sec,w3var,&         
                                        Skew_w,w1_1,w1_2,w2_1,w2_2,a)

  end subroutine shoc_assumed_pdf_vv_parameters_c   
  
  subroutine shoc_assumed_pdf_thl_parameters_c(&
                           wthlsec,sqrtw2,sqrtthl,thlsec,thl_first,& 
                           w1_1,w1_2,Skew_w,a,dothetal_skew,&
                           thl1_1,thl1_2,thl2_1,thl2_2,sqrtthl2_1,&
                           sqrtthl2_2,Skew_thl) bind (C)
    use shoc, only: shoc_assumed_pdf_thl_parameters

    real(kind=c_real), intent(in), value :: wthlsec
    real(kind=c_real), intent(in), value :: sqrtw2
    real(kind=c_real), intent(in), value :: sqrtthl
    real(kind=c_real), intent(in), value :: thlsec
    real(kind=c_real), intent(in), value :: thl_first
    real(kind=c_real), intent(in), value :: w1_1
    real(kind=c_real), intent(in), value :: w1_2
    real(kind=c_real), intent(in), value :: Skew_w
    real(kind=c_real), intent(in), value :: a
    logical(kind=c_bool), intent(in), value :: dothetal_skew

    real(kind=c_real), intent(out) :: thl1_1
    real(kind=c_real), intent(out) :: thl1_2
    real(kind=c_real), intent(out) :: thl2_1
    real(kind=c_real), intent(out) :: thl2_2
    real(kind=c_real), intent(out) :: sqrtthl2_1
    real(kind=c_real), intent(out) :: sqrtthl2_2
    real(kind=c_real), intent(out) :: Skew_thl

    call shoc_assumed_pdf_thl_parameters(&
                           wthlsec,sqrtw2,sqrtthl,thlsec,thl_first,& 
                           w1_1,w1_2,Skew_w,a,dothetal_skew,&
                           thl1_1,thl1_2,thl2_1,thl2_2,sqrtthl2_1,&
                           sqrtthl2_2,Skew_thl)

  end subroutine shoc_assumed_pdf_thl_parameters_c  
  
  subroutine shoc_assumed_pdf_qw_parameters_c(&
                           wqwsec,sqrtw2,Skew_w,sqrtqt,qwsec,& 
                           w1_1,w1_2,qw_first,a,&
                           qw1_1,qw1_2,qw2_1,qw2_2,sqrtqw2_1,&
                           sqrtqw2_2,Skew_qw) bind (C)
    use shoc, only: shoc_assumed_pdf_qw_parameters

    real(kind=c_real), intent(in), value :: wqwsec
    real(kind=c_real), intent(in), value :: sqrtw2
    real(kind=c_real), intent(in), value :: sqrtqt
    real(kind=c_real), intent(in), value :: qwsec
    real(kind=c_real), intent(in), value :: qw_first
    real(kind=c_real), intent(in), value :: w1_1
    real(kind=c_real), intent(in), value :: w1_2
    real(kind=c_real), intent(in), value :: Skew_w
    real(kind=c_real), intent(in), value :: a

    real(kind=c_real), intent(out) :: qw1_1
    real(kind=c_real), intent(out) :: qw1_2
    real(kind=c_real), intent(out) :: qw2_1
    real(kind=c_real), intent(out) :: qw2_2
    real(kind=c_real), intent(out) :: sqrtqw2_1
    real(kind=c_real), intent(out) :: sqrtqw2_2
    real(kind=c_real), intent(out) :: Skew_qw

    call shoc_assumed_pdf_qw_parameters(&
                          wqwsec,sqrtw2,Skew_w,sqrtqt,qwsec,& 
                          w1_1,w1_2,qw_first,a,&
                          qw1_1,qw1_2,qw2_1,qw2_2,sqrtqw2_1,&
                          sqrtqw2_2,Skew_qw)

  end subroutine shoc_assumed_pdf_qw_parameters_c
  
  subroutine shoc_assumed_pdf_inplume_correlations_c(&
                sqrtqw2_1,sqrtthl2_1,a,sqrtqw2_2,sqrtthl2_2,&
                qwthlsec,qw1_1,qw_first,thl1_1,thl_first,qw1_2,thl1_2,&
                r_qwthl_1) bind (C)
    use shoc, only: shoc_assumed_pdf_inplume_correlations

    real(kind=c_real), intent(in), value :: sqrtqw2_1
    real(kind=c_real), intent(in), value :: sqrtthl2_1
    real(kind=c_real), intent(in), value :: a
    real(kind=c_real), intent(in), value :: sqrtqw2_2
    real(kind=c_real), intent(in), value :: sqrtthl2_2
    real(kind=c_real), intent(in), value :: qwthlsec
    real(kind=c_real), intent(in), value :: qw1_1
    real(kind=c_real), intent(in), value :: qw_first
    real(kind=c_real), intent(in), value :: thl1_1
    real(kind=c_real), intent(in), value :: thl_first
    real(kind=c_real), intent(in), value :: qw1_2
    real(kind=c_real), intent(in), value :: thl1_2

    real(kind=c_real), intent(out) :: r_qwthl_1

    call shoc_assumed_pdf_inplume_correlations(&
                sqrtqw2_1,sqrtthl2_1,a,sqrtqw2_2,sqrtthl2_2,&
                qwthlsec,qw1_1,qw_first,thl1_1,thl_first,qw1_2,thl1_2,&
                r_qwthl_1)

  end subroutine shoc_assumed_pdf_inplume_correlations_c       

  subroutine shoc_pblintd_init_pot_c(shcol, nlev, thl, ql, q, thv) bind(C)
   use shoc, only: pblintd_init_pot

   integer(kind=c_int), value, intent(in) :: shcol, nlev
   real(kind=c_real), intent(in)  :: thl(shcol, nlev), ql(shcol, nlev), q(shcol, nlev)
   real(kind=c_real), intent(out) :: thv(shcol, nlev)

   call pblintd_init_pot(shcol, nlev, thl, ql, q, thv)

 end subroutine shoc_pblintd_init_pot_c


end module shoc_iso_c
