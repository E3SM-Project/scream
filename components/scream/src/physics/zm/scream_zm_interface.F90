#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module scream_zm_interface_mod

  use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_double, c_bool,C_NULL_CHAR, c_float
  use physics_utils, only: r8 => rtype, rtype8, itype, btype
  use zm_conv,       only: zm_convr
 
  implicit none
#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
  integer,parameter,public :: rtype = c_double ! 8 byte real, compatible with c type double
#else
# define c_real c_float
  integer,parameter,public :: rtype = c_float ! 4 byte real, compatible with c type #endif
#endif



  public :: zm_init_f90
  public :: zm_main_f90
  public :: zm_finalize_f90

  real   :: test
  
  real(r8) :: pcols = 32
  real(r8) :: pver = 72
  real(kind=c_real) :: cpair  !=    1004.64000000000
  real(kind=c_real) :: gravit !=    9.80616000000000
  real(kind=c_real) :: latvap !=    2501000.00000000
  real(kind=c_real) :: plevp = 0.0 

  real(kind=c_real) :: pverp = 0
  real(kind=c_real) :: ncnst = 0

contains

  !====================================================================!

  subroutine zm_init_f90 (pref_edge) bind(c)
    
    implicit none
    real(kind=c_real) :: pref_edge(plevp)  

    logical :: no_deep_pbl


    Print *, 'In zm_init_f90'

  end subroutine zm_init_f90
  !====================================================================!
  subroutine zm_main_f90 (lchnk   ,ncol    , &
                    t       ,qh      ,prec    ,jctop   ,jcbot   , &
                    pblh    ,zm      ,geos    ,zi      ,qtnd    , &
                    heat    ,pap     ,paph    ,dpp     , &
                    delt    ,mcon    ,cme     ,cape    , &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
                    mu      ,md      ,du      ,eu      ,ed      , &
                    dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
                    lengath ,ql      ,rliq    ,landfrac,hu_nm1  , &
                    cnv_nm1 ,tm1     ,qm1     ,t_star  ,q_star, dcape, q) bind(c)
   real(kind=c_real), intent(inout) :: t(pcols,pver) ! State array  kg/kg
   real(kind=c_real) :: qh(pcols,pver) 
   !
   ! real(kind-c_real), intent(in) :: zm(pcols,pver) 
   integer :: lchnk
   integer :: i,ii,j,k
   
   integer, intent(in) :: ncol
   real(r8), intent(out) :: prec(pcols)
   real(r8), intent(out) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8), intent(in) :: pblh(pcols)
   real(r8), intent(in) :: geos(pcols)
   real(r8), intent(in) :: zi(pcols,pver+1)
   real(r8), intent(out) :: qtnd(pcols,pver)           ! specific humidity tendency (kg/kg/s)
   real(r8), intent(out) :: heat(pcols,pver)           ! heating rate (dry static energy tendency, W/kg)
   real(r8), intent(out) :: mcon(pcols,pverp)
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)
   real(r8), intent(out) :: cape(pcols)        ! w  convective available potential energy.
   real(r8), intent(out) :: zdu(pcols,pver)
   real(r8), intent(out) :: rprd(pcols,pver)     ! rain production rate
! move these vars from local storage to output so that convective
! transports can be done in outside of conv_cam.
   real(r8), intent(out) :: mu(pcols,pver)
   real(r8), intent(in) :: zm(pcols,pver)
   real(r8), intent(in) :: pap(pcols,pver)     
   real(r8), intent(in) :: paph(pcols,pver+1)
   real(r8), intent(in) :: dpp(pcols,pver)        ! local sigma half-level thickness (i.e. dshj).
   real(r8), intent(in) :: tpert(pcols)
   real(r8), intent(out) :: md(pcols,pver)
   real(r8), intent(out) :: du(pcols,pver)       ! detrainement rate of updraft
   real(r8), intent(out) :: ed(pcols,pver)       ! entrainment rate of downdraft
   real(r8), intent(out) :: eu(pcols,pver)       ! entrainment rate of updraft
   real(r8), intent(out) :: dp(pcols,pver)       ! wg layer thickness in mbs (between upper/lower interface).
   real(r8), intent(out) :: dsubcld(pcols)       ! wg layer thickness in mbs between lcl and maxi.
   integer jt(pcols)                          ! wg top  level index of deep cumulus convection.
   integer maxg(pcols)                        ! wg gathered values of maxi.
   integer ideep(pcols)                       ! w holds position of gathered points vs longitude index.
   integer lengath
!     diagnostic field used by chem/wetdep codes
   real(r8), intent(out) :: ql(pcols,pver)
   real(r8) landfracg(pcols)            ! wg grid slice of landfrac  
   real(r8), intent(inout) :: hu_nm1 (pcols,pver)
   real(r8), intent(inout) :: cnv_nm1 (pcols,pver)
   real(r8), intent(in) :: tm1(pcols,pver)       ! grid slice of temperature at mid-layer.
   real(r8), intent(in) :: qm1(pcols,pver)       ! grid slice of specific humidity.
   real(r8), intent(in) :: landfrac(pcols) ! RBN Landfrac


   real(r8), intent(in), pointer, dimension(:,:) :: t_star ! intermediate T between n and n-1 time step
   real(r8), intent(in), pointer, dimension(:,:) :: q_star ! intermediate q between n and n-1 time step
   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
   real(r8), intent(out) :: dcape(pcols)           ! output dynamical CAPE

   real(r8) :: pcont(pcols), pconb(pcols), freqzm(pcols)
   real(r8), intent(out) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
   real(r8) delt                     ! length of model time-step in seconds.
   integer mx(pcols) 
   real(r8) q(pcols,pver)              ! w  grid slice of mixing ratio.
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   integer :: nstep             ! Time step index

   logical :: doconvtran(ncnst)
   real(r8) :: fake_dpdry(pcols,pver)       ! Delta pressure between interfaces


! input/output

   real(r8) :: fake_dqdt(pcols,pver,ncnst)  ! Tracer tendency array

   integer :: il1g
   fake_dpdry(:,:) = 0._r8

   fake_dqdt(:,:,1) = 0._r8
   il1g = 1
   doconvtran = .false.


   call zm_convr(lchnk   ,ncol, &
                    t       ,qh      ,prec    ,jctop   ,jcbot   , &
                    pblh    ,zm      ,geos    ,zi      ,qtnd    , &
                    heat    ,pap     ,paph    ,dpp     , &
                    delt    ,mcon    ,cme     ,cape    , &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
                    mu      ,md      ,du      ,eu      ,ed      , &
                    dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
                    lengath ,ql      ,rliq    ,landfrac,hu_nm1  , &
                    cnv_nm1 ,tm1     ,qm1     ,t_star  ,q_star, dcape) 
!   call convtran(lchnk   , &
!                    doconvtran,q       ,ncnst   ,mu      ,md      , &
!                    du      ,eu      ,ed      ,dp      ,dsubcld , &
!                    jt      ,mx      ,ideep   ,il1g    , lengath    , &
!                    1   ,fracis  ,fake_dqdt,  fake_dpdry   )

   end subroutine zm_main_f90
  !====================================================================!
  subroutine zm_finalize_f90 () bind(c)
    implicit none

  !  Print *, 'In zm_finalize_f90'

  end subroutine zm_finalize_f90
  !====================================================================!

end module scream_zm_interface_mod
