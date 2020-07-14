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
  integer,parameter,public :: rtype = c_float ! 4 byte real, compatible with c type float
#endif



  public :: zm_init_f90
  public :: zm_main_f90
  public :: zm_finalize_f90

  real   :: test
  
contains

  !====================================================================!
  subroutine zm_init_f90 () bind(c)
    
!    use zm_conv_intr,                   only: zm_conv_init
!    real(rtype), intent(in)  :: cpair  ! specific heat of dry air
!    call zm_conv_init(cpair)
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
