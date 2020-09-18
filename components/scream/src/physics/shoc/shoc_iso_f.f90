module shoc_iso_f
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from shoc fortran to scream c++.
!

interface

  subroutine calc_shoc_varorcovar_f(shcol, nlev, nlevi, tunefac, isotropy_zi, tkh_zi, dz_zi, invar1, invar2, varorcovar) bind (C)
    use iso_c_binding

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

  end subroutine calc_shoc_varorcovar_f

  subroutine calc_shoc_vertflux_f(shcol, nlev, nlevi, tkh_zi, dz_zi, invar, vertflux) bind (C)
    use iso_c_binding

    integer(kind=c_int), intent(in), value :: shcol
    integer(kind=c_int), intent(in), value :: nlev
    integer(kind=c_int), intent(in), value :: nlevi
    real(kind=c_real), intent(in) :: tkh_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: dz_zi(shcol,nlevi)
    real(kind=c_real), intent(in) :: invar(shcol,nlev)

    real(kind=c_real), intent(inout) :: vertflux(shcol,nlevi)

  end subroutine calc_shoc_vertflux_f

 subroutine shoc_diag_second_moments_srf_f(shcol, wthl, uw, vw, ustar2, wstar) bind(C)
   use iso_c_binding

   integer(kind=c_int), value, intent(in) :: shcol

   ! arguments
   real(kind=c_real), intent(in) :: wthl(shcol)
   real(kind=c_real), intent(in) :: uw(shcol)
   real(kind=c_real), intent(in) :: vw(shcol)
   real(kind=c_real), intent(out) :: ustar2(shcol)
   real(kind=c_real), intent(out) :: wstar(shcol)

 end subroutine shoc_diag_second_moments_srf_f

 subroutine shoc_diag_second_moments_ubycond_f(shcol, thl, qw, wthl, wqw, qwthl, uw, vw, wtke) bind(C)
   use iso_c_binding

   ! arguments
   integer(kind=c_int), value, intent(in) :: shcol
   real(kind=c_real), intent(out)  :: thl(shcol), qw(shcol), qwthl(shcol),wthl(shcol),wqw(shcol), &
        uw(shcol), vw(shcol), wtke(shcol)

 end subroutine shoc_diag_second_moments_ubycond_f

subroutine update_host_dse_f(shcol, nlev, thlm, shoc_ql, exner, zt_grid, &
                             phis, host_dse) bind (C)
  use iso_c_binding

  integer(kind=c_int), intent(in), value :: shcol
  integer(kind=c_int), intent(in), value :: nlev
  real(kind=c_real), intent(in) :: thlm(shcol,nlev)
  real(kind=c_real), intent(in) :: shoc_ql(shcol,nlev)
  real(kind=c_real), intent(in) :: exner(shcol,nlev)
  real(kind=c_real), intent(in) :: zt_grid(shcol,nlev)
  real(kind=c_real), intent(in) :: phis(shcol)

  real(kind=c_real), intent(out) :: host_dse(shcol,nlev)

end subroutine update_host_dse_f

 subroutine adv_sgs_tke_f(nlev, shcol, dtime, shoc_mix, wthv_sec, sterm_zt, tk, tke, a_diss) bind(C)
    use iso_c_binding

    !arguments
    integer(kind=c_int), value, intent(in) :: nlev, shcol

    real(kind=c_real),   value, intent(in)    :: dtime
    real(kind=c_real),          intent(in)    :: shoc_mix(shcol,nlev)
    real(kind=c_real),          intent(in)    :: wthv_sec(shcol,nlev)
    real(kind=c_real),          intent(in)    :: sterm_zt(shcol,nlev)
    real(kind=c_real),          intent(in)    :: tk(shcol,nlev)
    real(kind=c_real),          intent(inout) :: tke(shcol,nlev)
    real(kind=c_real),          intent(out)   :: a_diss(shcol,nlev)

 end subroutine adv_sgs_tke_f


end interface

end module shoc_iso_f
