module micro_p3_iso_f
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from micro_p3 fortran to scream c++.
!

interface

  subroutine find_lookuptable_indices_1a_f(dumi,dumjj,dumii,dumzz,dum1,dum4,dum5,dum6,      &
       qitot,nitot,qirim,rhop) bind(C)
    use iso_c_binding

    ! arguments:
    integer(kind=c_int), intent(out) :: dumi,dumjj,dumii,dumzz
    real(kind=c_real),   intent(out) :: dum1,dum4,dum5,dum6
    real(kind=c_real),   value, intent(in)  :: qitot,nitot,qirim,rhop
  end subroutine find_lookuptable_indices_1a_f

  subroutine find_lookuptable_indices_1b_f(dumj,dum3,qr,nr) bind(C)
    use iso_c_binding

    integer(kind=c_int), intent(out) :: dumj
    real(kind=c_real),   intent(out) :: dum3
    real(kind=c_real),   value, intent(in) :: qr, nr
  end subroutine find_lookuptable_indices_1b_f

  subroutine access_lookup_table_f(dumjj,dumii,dumi,index,dum1,dum4,dum5,proc) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: dumjj, dumii, dumi, index
    real(kind=c_real),   value, intent(in) :: dum1, dum4, dum5
    real(kind=c_real),   intent(out) :: proc
  end subroutine access_lookup_table_f

  subroutine access_lookup_table_coll_f(dumjj,dumii,dumj,dumi,index,dum1,dum3,dum4,dum5,proc) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: dumjj,dumii,dumj,dumi,index
    real(kind=c_real),   value, intent(in) :: dum1,dum3,dum4,dum5
    real(kind=c_real),   intent(out) :: proc
  end subroutine access_lookup_table_coll_f

  subroutine get_cloud_dsd2_f(qc,nc,mu_c,rho,nu,lamc,cdist,cdist1,lcldm) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)         :: qc,rho,lcldm
    real(kind=c_real), intent(inout)             :: nc
    real(kind=c_real), intent(out)               :: mu_c,nu,lamc,cdist,cdist1
  end subroutine get_cloud_dsd2_f

  subroutine get_rain_dsd2_f(qr,nr,mu_r,lamr,cdistr,logn0r,rcldm) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: qr,rcldm
    real(kind=c_real), intent(inout)     :: nr
    real(kind=c_real), intent(out)       :: lamr,mu_r,cdistr,logn0r
  end subroutine get_rain_dsd2_f

  subroutine cloud_water_autoconversion_f(rho, qc_incld, nc_incld, qcaut, ncautc, ncautr) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: rho, qc_incld, nc_incld
    real(kind=c_real), intent(inout) :: qcaut, ncautc, ncautr
  end subroutine cloud_water_autoconversion_f

  subroutine calc_first_order_upwind_step_f(kts, kte, kdir, kbot, k_qxtop, dt_sub, rho, inv_rho, inv_dzq, num_arrays, fluxes, vs, qnx) bind(C)
    use iso_c_binding

    !arguments:
    integer(kind=c_int), value, intent(in) :: kts, kte, kdir, kbot, k_qxtop, num_arrays
    real(kind=c_real), value, intent(in) :: dt_sub
    real(kind=c_real), dimension(kts:kte), intent(in) :: rho, inv_rho, inv_dzq
    type(c_ptr), intent(in), dimension(num_arrays) :: fluxes, vs, qnx
  end subroutine calc_first_order_upwind_step_f

  subroutine generalized_sedimentation_f(kts, kte, kdir, k_qxtop, k_qxbot, kbot, Co_max, dt_left, prt_accum, inv_dzq, inv_rho, rho, num_arrays, vs, fluxes, qnx) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: kts, kte, kdir, k_qxtop, kbot, num_arrays
    integer(kind=c_int), intent(inout) :: k_qxbot
    real(kind=c_real), value, intent(in) :: Co_max
    real(kind=c_real), intent(inout) :: dt_left, prt_accum
    real(kind=c_real), dimension(kts:kte), intent(in) :: inv_dzq, inv_rho, rho

    type(c_ptr), intent(in), dimension(num_arrays) :: vs, fluxes, qnx
  end subroutine generalized_sedimentation_f

  subroutine cloud_sedimentation_f(kts,kte,ktop,kbot,kdir,   &
       qc_incld,rho,inv_rho,lcldm,acn,inv_dzq,&
       dt,odt,log_predictNc, &
       qc, nc, nc_incld,mu_c,lamc,prt_liq,qc_tend,nc_tend) bind(C)

    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir

    real(kind=c_real), intent(in), dimension(kts:kte) :: qc_incld
    real(kind=c_real), intent(in), dimension(kts:kte) :: rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: lcldm
    real(kind=c_real), intent(in), dimension(kts:kte) :: acn
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_dzq

    real(kind=c_real),    value, intent(in) :: dt
    real(kind=c_real),    value, intent(in) :: odt
    logical(kind=c_bool), value, intent(in) :: log_predictNc

    real(kind=c_real), intent(inout), dimension(kts:kte) :: qc
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc_incld
    real(kind=c_real), intent(inout), dimension(kts:kte) :: mu_c
    real(kind=c_real), intent(inout), dimension(kts:kte) :: lamc
    real(kind=c_real), intent(inout) :: prt_liq
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qc_tend
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nc_tend
  end subroutine cloud_sedimentation_f

  subroutine ice_sedimentation_f(kts,kte,ktop,kbot,kdir,    &
       rho,inv_rho,rhofaci,icldm,inv_dzq,dt,odt,  &
       qitot,qitot_incld,nitot,qirim,qirim_incld,birim,birim_incld,nitot_incld,prt_sol,qi_tend,ni_tend) bind(C)

    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir

    real(kind=c_real), intent(in), dimension(kts:kte) :: rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: rhofaci
    real(kind=c_real), intent(in), dimension(kts:kte) :: icldm
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_dzq
    real(kind=c_real), value, intent(in) :: dt, odt

    real(kind=c_real), intent(inout), dimension(kts:kte), target :: qitot
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qitot_incld
    real(kind=c_real), intent(inout), dimension(kts:kte), target :: nitot
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nitot_incld
    real(kind=c_real), intent(inout), dimension(kts:kte), target :: qirim
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qirim_incld
    real(kind=c_real), intent(inout), dimension(kts:kte), target :: birim
    real(kind=c_real), intent(inout), dimension(kts:kte) :: birim_incld

    real(kind=c_real), intent(inout) :: prt_sol
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qi_tend
    real(kind=c_real), intent(inout), dimension(kts:kte) :: ni_tend
  end subroutine ice_sedimentation_f

  subroutine rain_sedimentation_f(kts,kte,ktop,kbot,kdir,   &
       qr_incld,rho,inv_rho,rhofacr,rcldm,inv_dzq,dt,odt,  &
       qr,nr,nr_incld,mu_r,lamr,prt_liq,rflx,qr_tend,nr_tend) bind(C)
    use iso_c_binding

    integer(kind=c_int), value, intent(in) :: kts, kte, ktop, kbot, kdir

    real(kind=c_real), intent(in), dimension(kts:kte) :: qr_incld

    real(kind=c_real), intent(in), dimension(kts:kte) :: rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_rho
    real(kind=c_real), intent(in), dimension(kts:kte) :: rhofacr
    real(kind=c_real), intent(in), dimension(kts:kte) :: rcldm
    real(kind=c_real), intent(in), dimension(kts:kte) :: inv_dzq
    real(kind=c_real), value, intent(in) :: dt, odt

    real(kind=c_real), intent(inout), target, dimension(kts:kte) :: qr
    real(kind=c_real), intent(inout), target, dimension(kts:kte) :: nr
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nr_incld
    real(kind=c_real), intent(inout), dimension(kts:kte) :: mu_r
    real(kind=c_real), intent(inout), dimension(kts:kte) :: lamr
    real(kind=c_real), intent(inout) :: prt_liq
    real(kind=c_real), intent(inout), dimension(kts:kte+1) :: rflx
    real(kind=c_real), intent(inout), dimension(kts:kte) :: qr_tend
    real(kind=c_real), intent(inout), dimension(kts:kte) :: nr_tend

  end subroutine rain_sedimentation_f

  subroutine calc_bulk_rho_rime_f(qi_tot, qi_rim, bi_rim, rho_rime) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real),   value, intent(in)  :: qi_tot
    real(kind=c_real),   intent(inout) :: qi_rim, bi_rim
    real(kind=c_real),   intent(out) :: rho_rime
  end subroutine calc_bulk_rho_rime_f

  subroutine compute_rain_fall_velocity_f(qr_incld, rcldm, rhofacr, nr, nr_incld, mu_r, lamr, V_qr, V_nr) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), value, intent(in) :: qr_incld, rcldm, rhofacr
    real(kind=c_real), intent(inout) :: nr, nr_incld
    real(kind=c_real), intent(out) :: mu_r, lamr, V_qr, V_nr
  end subroutine compute_rain_fall_velocity_f

  subroutine ice_cldliq_collection_f(rho, t, rhofaci, f1pr04, qitot_incld, qc_incld, nitot_incld, &
                                     nc_incld, qccol, nccol, qcshd, ncshdc) bind(C)
    use iso_c_binding
    
    ! arguments:
    real(kind=c_real), value, intent(in) :: rho, t, rhofaci, f1pr04
    real(kind=c_real), intent(in) :: qitot_incld, qc_incld, nitot_incld, nc_incld
    real(kind=c_real), intent(out) :: qccol, nccol, qcshd, ncshdc
  end subroutine ice_cldliq_collection_f

  subroutine ice_rain_collection_f(rho, t, rhofaci, logn0r, f1pr07, f1pr08, &
                                   qitot_incld, nitot_incld, qr_incld, qrcol, nrcol) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), intent(in) :: rho, t, rhofaci, logn0r, f1pr07, f1pr08
    real(kind=c_real), intent(in) :: qitot_incld, nitot_incld, qr_incld
    real(kind=c_real), intent(out) :: qrcol, nrcol
  end subroutine ice_rain_collection_f

  subroutine ice_self_collection_f(rho, rhofaci, f1pr03, eii, qirim_incld, qitot_incld, nitot_incld, nislf) bind(C)
    use iso_c_binding

    ! arguments:
    real(kind=c_real), intent(in) :: rho, rhofaci, f1pr03, eii, qirim_incld, qitot_incld, nitot_incld
    real(kind=c_real), intent(out) :: nislf
  end subroutine ice_self_collection_f 
  !
  ! These are some routine math operations that are not BFB between
  ! fortran and C++ on all platforms, so fortran will need to use
  ! the C++ versions in order to stay BFB.
  !

  function cxx_pow(base, exp) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base
    real(kind=c_real), value, intent(in)  :: exp

    ! return
    real(kind=c_real)               :: cxx_pow
  end function cxx_pow

  function cxx_cbrt(base) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base

    ! return
    real(kind=c_real)               :: cxx_cbrt
  end function cxx_cbrt

  function cxx_gamma(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_gamma
  end function cxx_gamma

  function cxx_log(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_log
  end function cxx_log

  function cxx_log10(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_log10
  end function cxx_log10

  function cxx_exp(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_exp
  end function cxx_exp

end interface

end module micro_p3_iso_f
