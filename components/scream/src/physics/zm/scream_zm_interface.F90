#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module scream_zm_interface_mod

  use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_double, c_bool,C_NULL_CHAR, c_float
 
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

  integer(kind=c_int) :: pcols = 32
  integer(kind=c_int) :: pver  = 72
!  integer(kind=c_int) :: ncol

  real(kind=c_real) :: cpair  !=    1004.64000000000
  real(kind=c_real) :: gravit !=    9.80616000000000
  real(kind=c_real) :: latvap !=    2501000.00000000
  
  real(kind=c_real) :: plevp = 0.0 

contains

  !====================================================================!
  subroutine zm_init_f90 () bind(c)
    use zm_conv,        only: zm_convi
 
    implicit none
    real(kind=c_real) :: pref_edge(plevp)  

    logical :: no_deep_pbl

    real(kind=c_real) :: prec !total precipitation from ZM convection
    real(kind=c_real) :: mdt !  tendency - Zhang-McFarlane moist convection
    real(kind=c_real) :: mdq ! Zhang-McFarlane moist convection
    real(kind=c_real) :: cice ! Cloud ice tendency - Zhang-McFarlane convection
    real(kind=c_real) :: cliq ! Cloud liq tendency - Zhang-McFarlane convection
    real(kind=c_real) :: evap ! Evaporation/snow prod from Zhang convection
    real(kind=c_real) :: rain_snow_convr ! Rain to snow conversion from Zhang convection
    real(kind=c_real) :: snow_rain_convr ! Snow to rain prod from Zhang convection
    real(kind=c_real) :: evapmoistcnv ! Evaporation from Zhang-McFarlane moist convection
    real(kind=c_real) :: pflux ! Flux of precipitation from ZM convection       
    real(kind=c_real) :: sflux ! Flux of snow from ZM convection                
    real(kind=c_real) :: netp ! Net precipitation production from ZM convection
    real(kind=c_real) :: netsn ! Net snow production from ZM convection         
    real(kind=c_real) :: heat ! Heating by ice and evaporation in ZM convection'
    real(kind=c_real) :: massflux ! Convection mass flux from ZM deep 
    real(kind=c_real) :: prrate ! Convective precipitation rate from ZM deep

    real(kind=c_real) :: basepres ! convection base pressure
    real(kind=c_real) :: tppres ! convection top  pressure

    real(kind=c_real) :: llvl ! model level of launching parcel

    real(kind=c_real) :: av_pe ! Convectively available potential energy
    real(kind=c_real) :: conv_occr ! Fractional occurance of ZM convection

   
    real(kind=c_real) :: tend ! T tendency - ZM convective momentum transport
    real(kind=c_real) :: tendt ! T tendency - ZM convective momentum transport
    real(kind=c_real) :: tendv ! V tendency - ZM convective momentum transport

    real(kind=c_real) :: upflx ! ZM convection updraft mass flux
    real(kind=c_real) :: downflx ! ZM convection downdraft mass flux'

    real(kind=c_real) :: znl_up_frc ! zonal force from ZM updraft pressure gradient term
    real(kind=c_real) :: znl_dwn_frc ! zonal force from ZM downdraft pressure gradient term
    real(kind=c_real) :: mrdl_up_frc ! meridional force from ZM updraft pressure gradient term
    real(kind=c_real) :: mrdl_dwn_frc ! merdional force from ZM downdraft pressure gradient term
    real(kind=c_real) :: uupd ! ZM in-cloud U updrafts
    real(kind=c_real) :: udwnd  ! ZM in-cloud U downdrafts
    real(kind=c_real) :: vupd ! ZM in-cloud V updrafts
    real(kind=c_real) :: vdwnd ! ZM in-cloud V downdrafts
    Print *, 'In zm_init_f90'

  end subroutine zm_init_f90
  !====================================================================!
  subroutine zm_main_f90 () bind(c)
    implicit none
    
 !   Print *, 'In zm_main_f90 '

  end subroutine zm_main_f90
  !====================================================================!
  subroutine zm_finalize_f90 () bind(c)
    implicit none

  !  Print *, 'In zm_finalize_f90'

  end subroutine zm_finalize_f90
  !====================================================================!

end module scream_zm_interface_mod
