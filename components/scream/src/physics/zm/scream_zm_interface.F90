#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module scream_zm_interface_mod

  use iso_c_binding, only: c_ptr, c_f_pointer, c_int, c_double, c_bool,C_NULL_CHAR, c_float

  implicit none


  public :: zm_init_f90
  public :: zm_main_f90
  public :: zm_finalize_f90

  real   :: test
  

contains

  !====================================================================!
  subroutine zm_init_f90 () bind(c)

    Print *, 'In zm_init_f90'

  end subroutine zm_init_f90
  !====================================================================!
  subroutine zm_main_f90 () bind(c)
    
    Print *, 'In zm_main_f90 '

  end subroutine zm_main_f90
  !====================================================================!
  subroutine zm_finalize_f90 () bind(c)

    Print *, 'In zm_finalize_f90'

  end subroutine zm_finalize_f90
  !====================================================================!

end module scream_zm_interface_mod
