module read_spa_data

  !====================================================================================
  ! These set of routines uses tracer data codes to read and interpolate (space and
  ! time) spa (simplified prescribed aerosol) input data which is read in from
  ! an input file

  ! Author: Balwinder Singh
  !====================================================================================

  use shr_kind_mod, only: shr_kind_cl
  use tracer_data,  only: trfld, trfile ! data strutures for storing file and fields info
  use radconstants, only: nswbands, nlwbands ! short and long wave band numbers
  implicit none

  !declare all variables private by default
  private

  !list of public subroutines
  public :: read_spa_data_init, read_spa_data_adv, read_spa_data_register, spa_readnl

  !---------------------------------------------------------------------------
  !Module level data shared by more than one subroutine or protected variables
  !---------------------------------------------------------------------------

  logical, public, protected :: is_spa_active = .false.

  type(trfld), pointer :: spa_fields_type(:)
  type(trfile)         :: spa_file_type
  logical              :: rmv_file = .false.

  !Following fields are read in via namelist but defaults are set here
  character(len=shr_kind_cl) :: filename = ' '
  character(len=shr_kind_cl) :: datapath = ' '
  character(len=32)          :: datatype = 'SERIAL'
  integer                    :: cycle_yr = 0

  !fields to read from the file
character(len=16), parameter :: pbuf_names(59) = &
     [ 'AER_G_SW_0   ', 'AER_G_SW_1   ', 'AER_G_SW_2   ', 'AER_G_SW_3   ', 'AER_G_SW_4   ', &
       'AER_G_SW_5   ', 'AER_G_SW_6   ', 'AER_G_SW_7   ', 'AER_G_SW_8   ', 'AER_G_SW_9   ', &
       'AER_G_SW_10  ', 'AER_G_SW_11  ', 'AER_G_SW_12  ', 'AER_G_SW_13  ', &
       'AER_SSA_SW_0 ', 'AER_SSA_SW_1 ', 'AER_SSA_SW_2 ', 'AER_SSA_SW_3 ', 'AER_SSA_SW_4 ', &
       'AER_SSA_SW_5 ', 'AER_SSA_SW_6 ', 'AER_SSA_SW_7 ', 'AER_SSA_SW_8 ', 'AER_SSA_SW_9 ', &
       'AER_SSA_SW_10', 'AER_SSA_SW_11', 'AER_SSA_SW_12', 'AER_SSA_SW_13', &
       'AER_TAU_LW_0 ', 'AER_TAU_LW_1 ', 'AER_TAU_LW_2 ', 'AER_TAU_LW_3 ', 'AER_TAU_LW_4 ', &
       'AER_TAU_LW_5 ', 'AER_TAU_LW_6 ', 'AER_TAU_LW_7 ', 'AER_TAU_LW_8 ', 'AER_TAU_LW_9 ', &
       'AER_TAU_LW_10', 'AER_TAU_LW_11', 'AER_TAU_LW_12', 'AER_TAU_LW_13', 'AER_TAU_LW_14', &
       'AER_TAU_LW_15', &
       'AER_TAU_SW_0 ', 'AER_TAU_SW_1 ', 'AER_TAU_SW_2 ', 'AER_TAU_SW_3 ', 'AER_TAU_SW_4 ', &
       'AER_TAU_SW_5 ', 'AER_TAU_SW_6 ', 'AER_TAU_SW_7 ', 'AER_TAU_SW_8 ', 'AER_TAU_SW_9 ', &
       'AER_TAU_SW_10', 'AER_TAU_SW_11', 'AER_TAU_SW_12','AER_TAU_SW_13', &
       'CCN3         ']

character(len=16), parameter :: specifier(59) = pbuf_names(59)


contains

  !-------------------------------------------------------------------
  ! registers spa fields to the phys buffer
  !-------------------------------------------------------------------
  subroutine read_spa_data_register()

    use ppgrid,         only: pver,pcols
    use physics_buffer, only : pbuf_add_field, dtype_r8

    implicit none

    integer :: i,idx

    do i = 1,size(pbuf_names)
       call pbuf_add_field(pbuf_names(i),'physpkg',dtype_r8,(/pcols,pver/),idx)
    enddo

  endsubroutine read_spa_data_register

  !-------------------------------------------------------------------
  ! Read spa namelist
  !-------------------------------------------------------------------
  subroutine spa_readnl(nlfile)

    use shr_log_mod,     only: errMsg => shr_log_errMsg
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use spmd_utils,      only: masterproc
    use mpishorthand,    only: mpiint, mpichar, mpicom
    use cam_abortutils,   only : endrun

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr

    character(len=shr_kind_cl) :: spa_file
    character(len=shr_kind_cl) :: spa_datapath
    character(len=32)          :: spa_type
    integer                    :: spa_cycle_yr

    !namelist fields to read
    namelist /spa_nl/ &
         spa_file,      & !file containing spa data
         spa_datapath,  & !path to the file
         spa_type,      & !Type of data(CYCLICAL,SERIAL,INTERP_MISSING_MONTHS,FIXED)
         spa_cycle_yr     !Year to cycle for the data

    !Initialize namelist variables from local module variables.
    spa_file     = filename
    spa_datapath = datapath
    spa_type     = datatype
    spa_cycle_yr = cycle_yr !BALLI replace cyle year with some other var name

    !Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'spa_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, spa_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('ERROR reading spa namelist, '//errmsg(__FILE__,__LINE__))
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
   !Broadcast namelist variables
   call mpibcast(spa_file,     len(spa_file),     mpichar, 0, mpicom)
   call mpibcast(spa_datapath, len(spa_datapath), mpichar, 0, mpicom)
   call mpibcast(spa_type,     len(spa_type),     mpichar, 0, mpicom)
   call mpibcast(spa_cycle_yr, 1, mpiint,  0, mpicom)
#endif

   !Update module variables with user settings.
   filename   = spa_file
   datapath   = spa_datapath
   datatype   = spa_type
   cycle_yr   = spa_cycle_yr

   !Turn on spa if user has specified an input dataset.
   is_spa_active = len_trim(filename) > 0

  end subroutine spa_readnl

  !-------------------------------------------------------------------
  ! Initialize spa input data to be read from an input file
  !-------------------------------------------------------------------
  subroutine read_spa_data_init

    use tracer_data, only: trcdata_init

    implicit none

    !Tracer data routine init
    allocate (spa_file_type%in_pbuf(size(specifier)))
    spa_file_type%in_pbuf(:) = .true.

!    call trcdata_init( specifier, 'unfied_SPA_file_lat_lon.nc', '', '/compyfs/sing201/lat_lon', spa_fields_type, spa_file_type, &
!         rmv_file, 1, 0, 0, 'CYCLICAL')
    call trcdata_init( specifier, 'unfied_SPA_file_lat_lon.nc', '', '/compyfs/sing201/lat_lon', spa_fields_type, spa_file_type, &
         rmv_file, 1, 0, 0, 'CYCLICAL')


  end subroutine read_spa_data_init

  !-------------------------------------------------------------------
  ! Advance fields in time and interpolate both in space and time
  !-------------------------------------------------------------------
  subroutine read_spa_data_adv( state, pbuf2d )

    !advance fields in time and interpolate (space and time)
    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk

    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !if(.not. is_spa_active) return

    !interpolate in time and space
    call advance_trcdata( spa_fields_type, spa_file_type, state, pbuf2d )

  end subroutine read_spa_data_adv

end module read_spa_data
