module nucleate_ice_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for nucleate_ice module.
!
!  B. Eaton - Sept 2014
!---------------------------------------------------------------------------------


use ppgrid,         only: pcols, pver
use physics_buffer, only: pbuf_add_field, dtype_r8


implicit none
private
public :: nucleate_ice_cam_register

integer :: naai_idx

!===============================================================================
contains
!===============================================================================

subroutine nucleate_ice_cam_register()



end subroutine nucleate_ice_cam_register

end module nucleate_ice_cam
