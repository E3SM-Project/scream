
module microp_aero

use ppgrid,           only: pver

implicit none
private
save

public ::  microp_aero_register
integer :: npccn_idx
integer :: naai_idx

contains
subroutine microp_aero_register
   use ppgrid,         only: pcols
   use physics_buffer, only: pbuf_add_field, dtype_r8

   call pbuf_add_field('NPCCN',      'physpkg',dtype_r8,(/pcols,pver/), npccn_idx)
   call pbuf_add_field('NAAI',     'physpkg', dtype_r8, (/pcols,pver/), naai_idx)
end subroutine microp_aero_register

end module microp_aero
