module nvtx

  implicit none

  !character(len=256),private :: tempName
  !character(len=256),private,target :: tempName
  !CHARACTER(c_char), ALLOCATABLE, TARGET :: tempName(:)

  interface nvtxRangePush
     subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
       use iso_c_binding, only : c_char
       character(kind=c_char,len=*) :: name
     end subroutine nvtxRangePushA
  end interface nvtxRangePush

  interface nvtxRangePop
     subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
     end subroutine nvtxRangePop
  end interface nvtxRangePop

contains

  subroutine nvtxStartRange(name,id)
    use iso_c_binding, only : c_char, c_loc, c_null_char
    !character(kind=c_char,len=*), intent(in) :: name
    character(len=*), intent(in) :: name
    integer, optional:: id
    !character(len=256) :: tempName
    !tempName=trim(name)//c_null_char

    !write(*,*) "ndk nvtxStartRange trim(name)=", trim(name)
    !call nvtxRangePush(tempName)
    !call nvtxRangePush(name)
!    call nvtxRangePush(trim(name)//c_null_char)
    !call nvtxRangePush(trim(name))

  end subroutine nvtxStartRange

  subroutine nvtxEndRange
!    call nvtxRangePop
  end subroutine nvtxEndRange

end module nvtx
