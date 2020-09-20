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

 subroutine shoc_pblintd_init_pot_f(shcol, nlev, thl, ql, q, thv) bind(C)
   use iso_c_binding

   integer(kind=c_int), value, intent(in) :: shcol, nlev
   real(kind=c_real), intent(in)  :: thl(shcol, nlev), ql(shcol, nlev), q(shcol,nlev)
   real(kind=c_real), intent(out) :: thv(shcol, nlev)

 end subroutine shoc_pblintd_init_pot_f

 subroutine shoc_diag_second_moments_ubycond_f(shcol, thl, qw, wthl, wqw, qwthl, uw, vw, wtke) bind(C)
   use iso_c_binding

   ! argmens
   integer(kind=c_int), value, intent(in) :: shcol
   real(kind=c_real), intent(out)  :: thl(shcol), qw(shcol), qwthl(shcol),wthl(shcol),wqw(shcol), &
        uw(shcol), vw(shcol), wtke(shcol)

 end subroutine shoc_diag_second_moments_ubycond_f

end interface

end module shoc_iso_f
