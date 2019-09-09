#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module simple_io_mod
!=====================================================================!
! Fortran/C++ interface for handling simple input and output.
! 
! Developed in order to allow us to load input files and write output
! to be used for convergence and baseline testing.
!
! The plan is to DELETE this interface when the more sophisticated
! PIO2 interface is adopted by SCREAM.
!
! TODO: Replace this I/O interface with PIO2
!
! Author: Aaron S. Donahue
! email: donahue5@llnl.gov
! data: Septemper 9, 2019
!=====================================================================!

  use netcdf
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

  logical(kind=c_bool) :: isopen_out = .false.
  logical(kind=c_bool) :: isopen_in  = .false.
  integer(kind=c_int)  :: max_ndims = 0

  integer(kind=c_int), public :: ncid_out  ! netcdf file id for output
  integer(kind=c_int), public :: ncid_in  ! netcdf file id for input
  integer(kind=c_int), allocatable :: start(:)

contains

!======================================================================
! TODO: Delete this, it is only for testing purposes at the moment.
  subroutine hello_world_f90(myint,myreal,mychar_c,len_planet,planets_c)  bind(c)
    integer(kind=c_int), value, intent(in) :: myint
    real(kind=c_real), value, intent(in) :: myreal
    type(c_ptr), intent(in) :: mychar_c
    integer(kind=c_int), value, intent(in) :: len_planet
    type(c_ptr), intent(in) :: planets_c(len_planet)
    character(kind=c_char, len=128) :: string
    
    integer :: len, lenp, i
    character(len=256) :: mychar
    character(len=256) :: planets
    character(len=256), dimension(len_planet) :: myplanets

    call convert_c_string(mychar_c,mychar)
    call convert_c_strarray(len_planet,planets_c,myplanets)

    string = "Hello World"
    print *, ""
    print *, trim(string), myint, myreal, trim(mychar)
    do i = 1,len_planet
      print *, 'Hello ', trim(myplanets(i))
    end do
    print *, ""

  end subroutine hello_world_f90
!======================================================================
  subroutine simple_in_init(filename_in) bind(c)
  ! Purpose: Open a netcdf file to be used for input to C++.
  ! Recieves a c-string filename which is converted to a fortran string
  ! and passed to the netcdf file open subroutine.

    type(c_ptr), intent(in) :: filename_in ! INPUT: Name of input file.
    
    character(len=256)       :: filename   ! Local Fortran variable to store filename in proper fortran string format.

    call convert_c_string(filename_in,filename) ! Convert c-string to fortran-string
    call check( nf90_open(filename, nf90_nowrite, ncid_in), "Open input file" ) ! Open netCDF file for input
    isopen_in = .true. ! Register that an input file has been opened.

    return
  end subroutine simple_in_init
!======================================================================
  subroutine simple_in_finalize() bind(c)
  ! Purpose: Close a netcdf file  used for input.

    call check( nf90_close(ncid_in), "Close input file" ) ! Close the input netcdf file
    isopen_in = .false. ! Register that the input file is now closed.

    return
  end subroutine simple_in_finalize
!======================================================================
  subroutine simple_io_get_field_size(fieldname_in,ndims,ldims) bind(c)

    type(c_ptr), intent(in)          :: fieldname_in
    integer(kind=c_int), value, intent(in)  :: ndims
    integer(kind=c_int), intent(out) :: ldims(ndims)

    character(len=256) :: fieldname
    integer :: varid, ndim_local, dimids(ndims)
    integer :: n
  
    call convert_c_string(fieldname_in,fieldname)
    
    call check( nf90_inq_varid(ncid_in, fieldname, varid), "get VAR Id" )
    call check( nf90_inquire_variable(ncid_in, varid, dimids=dimids, ndims=ndim_local), "get VAR-Dim Ids" )

    ldims(:) = 1
    do n = 1,ndim_local
      call check( nf90_inquire_dimension(ncid_in, dimids(n), len=ldims(n)), "get Var-Dim Len's" )
    end do
    
    return
  end subroutine simple_io_get_field_size
!======================================================================
  subroutine simple_out_init1(filename_in,ndims,dimnames_in,dimrng) bind(c)
  ! Purpose: To open (and replace if needed) a netCDF file to be used for model
  ! output.
  ! Receives a c-string filename which is used to create the new file.
  ! ndims is the total number of potential dimensions we anticpate needing for the output file 
  ! dimnames_in is a C++ str::array that contains the names for each dimension (len=ndims)
  ! dimrng is a C++ integer vector with the len of each dimension.
  
    type(c_ptr), intent(in) :: filename_in                    ! Name of netCDF output file
    integer(kind=c_int), value, intent(in)  :: ndims          ! Number of potential dimensions
    type(c_ptr), intent(in) :: dimnames_in(ndims)             ! Name for each dimension
    integer(kind=c_int), target, intent(in)  :: dimrng(ndims) ! Max range for each dimension, 0 for unlimited
  
    character(len=256) :: dimnames(ndims) ! Local string array of dimension names in fortran format
    character(len=256) :: filename        ! Local string of output filename in fortran format
  
    integer(kind=c_int) :: dimids(ndims)  ! The netcdf id for each dimension
    
    integer(kind=c_int) :: ii
 
    ! Convert C++ strings to fortran format 
    call convert_c_string(filename_in,filename)
    call convert_c_strarray(ndims,dimnames_in,dimnames)
  
    ! Open the netCDF output file (replace if needed)
    call check( nf90_create(trim(filename), nf90_clobber, ncid_out), 'init1' )
    isopen_out = .true.  ! Register that an output file is open
  
    ! Add each dimension to the netCDF output.
    do ii = 1,ndims
      call check( nf90_def_dim(ncid_out, trim(dimnames(ii)), dimrng(ii), dimids(ii)),'add dimensions' )
    end do
  
    return
  end subroutine simple_out_init1
!======================================================================
  subroutine simple_out_init2() bind(c)
  ! Purpose: To finish the registration of output variables.
  
    allocate(start(max_ndims))  ! TODO: Come up with a better way to do this. Possible return the max_ndims as an output
  
    ! Finish the definition phase of the netCDF output file
    call check( nf90_enddef(ncid_out), 'init2' )
    isopen_out = .false.
  
  end subroutine simple_out_init2
!======================================================================
  subroutine simple_out_finalize() bind(c)
  ! Purpose: To close the netCDF output file.
  
    call check( nf90_close(ncid_out), 'finalize' )
  
  end subroutine simple_out_finalize
!======================================================================
  subroutine simple_out_regfield(field_name_in,field_type,ndim,field_dim_in,units_in) bind(c)
  ! Purpose: To register an output field with the netCDF output file.
  ! Receives a C++ string, field_name_in, which will be the field label.
  ! Needs a field_type designation which is an integer id for variable type
  ! (such as real, int, char)
  ! ndim is an integer specifying the dimensionality of the field.
  ! field_dim_in is a str::array of dimensions for the field (which should match
  ! one of the dimensions declared in the first init phase.
  ! units_in is an optional C++ string that will declare the units for the
  ! field.
  
    type(c_ptr),         intent(in)           :: field_name_in       ! The field name
    integer(kind=c_int), value, intent(in)    :: field_type          ! The variable type (i.e. real, int, char)
    integer(kind=c_int), value, intent(in)    :: ndim                ! The dimensionality of the field (i.e. 1, 2, 3)
    type(c_ptr),         intent(in)           :: field_dim_in(ndim)  ! The dimensions of the field
    type(c_ptr),         intent(in), optional :: units_in            ! The units for the field (optional)
  
    character(len=256) :: field_name      ! Fortran string for field name
    character(len=256) :: field_dim(ndim) ! Fortran string array for dimension names
    character(len=256) :: units           ! Fortran sting for units
  
    integer(kind=c_int) :: varid          ! Local int to store the variable id from netCDF output
    integer(kind=c_int) :: dimids(ndim)   ! Local array to store dimension ids
    character(kind=c_char,len=100) :: check_case ! Error message 
    integer(kind=c_int) :: ii
  
    ! Make sure file is still open for variable registration
    if (.not.isopen_out) then
      print *, 'netCDF definition is closed, cannot register new field', field_name
      stop
    end if
  
    ! Convert C++ strings to fortran format
    call convert_c_string(field_name_in,field_name)
    call convert_c_strarray(ndim,field_dim_in,field_dim)
    if (present(units_in)) call convert_c_string(units_in,units)
  
    ! Get dimension ID's for dimensions related to field
    do ii = 1,ndim
      call check( nf90_inq_dimid(ncid_out,trim(field_dim(ii)),dimids(ii)), 'Reg: get dim ids' )
    end do
  
    ! Register the variable and units if applicable.
    check_case = 'Add variable '//trim(field_name)
    call check( nf90_def_var(ncid_out,trim(field_name), field_type, dimids, varid), check_case )
    check_case = 'Add variable '//trim(field_name)//', units'
    if (present(units_in)) then
      call check( nf90_put_att(ncid_out, varid, "units", trim(units)), check_case )
    else
      call check( nf90_put_att(ncid_out, varid, "units", "unspecified"), check_case )
    end if 
  
    ! Update the max number of dimensions for output file.
    max_ndims = max(max_ndims,ndim)
  
    return  
  end subroutine simple_out_regfield
!======================================================================
   subroutine io_readfield_real(field_name_in, time_dim, flen, field_data) bind(c)
   
   type(c_ptr), intent(in)      :: field_name_in
   integer(kind=c_int), intent(in) :: time_dim(2)
   integer(kind=c_int), value, intent(in) :: flen
   real(kind=c_real),intent(out) :: field_data(flen)

   real(kind=c_real), allocatable :: local_data(:,:,:,:,:,:)
   character(len=256)  :: field_name
   integer             :: ndims, varid
   integer             :: dimids(6), dimlens(6), start(6)
   integer(kind=c_int) :: i, dimprod

   call convert_c_string(field_name_in,field_name)

   call check( nf90_inq_varid(ncid_out, trim(field_name), varid), 'Get variable ID' )
   call check( nf90_inquire_variable(ncid_in,varid,ndims = ndims,dimids=dimids), 'Get variable DIMS' )
   
   dimprod = 1 ! Check that flen and dim lens match
   dimlens(:) = 0
   do i = 1,ndims
     call check( nf90_inquire_dimension(ncid_in, dimids(i), len=dimlens(i)),  'Get dim lens' ) 
     dimprod = dimprod*dimlens(i)
   end do
   if (dimprod.ne.flen) then
     print *, 'Field len = ', flen, ' and total dimension len = ',dimprod, ' does not match, exiting...'
     stop
   end if

   select case (ndims)
     case (1)
       allocate(local_data(dimlens(1),1,1,1,1,1))
     case (2)
       allocate(local_data(dimlens(1),dimlens(2),1,1,1,1))
     case (3)
       allocate(local_data(dimlens(1),dimlens(2),dimlens(3),1,1,1))
     case (4)
       allocate(local_data(dimlens(1),dimlens(2),dimlens(3),dimlens(4),1,1))
     case (5)
       allocate(local_data(dimlens(1),dimlens(2),dimlens(3),dimlens(4),dimlens(5),1))
     case (6)
       allocate(local_data(dimlens(1),dimlens(2),dimlens(3),dimlens(4),dimlens(5),dimlens(6)))
     case default
       print *, "READFIELD: ndims > 6 not supported"
       stop
   end select
   start(:) = 1
   start(time_dim(1)) = time_dim(2)
   call check( nf90_get_var(ncid_in,varid, local_data, start=start), 'Get variable info')
   field_data = reshape(local_data, (/ flen /))

   return
   end subroutine io_readfield_real
   !******************************************************************
   ! 1D
   subroutine io_readfield_1d_real(field_name_in,time_dim,dims,field_data) bind(c)

   type(c_ptr), intent(in)      :: field_name_in
   integer(kind=c_int), value, intent(in) :: time_dim
   integer(kind=c_int), value, intent(in) :: dims
   real(kind=c_real),intent(out) :: field_data(dims)

   character(len=256)  :: field_name
   integer(kind=c_int) :: ndims
   integer(kind=c_int) :: varid
   integer(kind=c_int) :: this_dim

   call convert_c_string(field_name_in,field_name)
   this_dim = 1 

   call check( nf90_inq_varid(ncid_out, trim(field_name), varid), 'Get variable ID' )
   call check( nf90_inquire_variable(ncid_out,varid,ndims = ndims), 'Get variable DIMS' )
   if (ndims.gt.this_dim) then ! Then there is a time dimension, handle appropriately.
     start(:) = 1
     start(ndims) = time_dim
     call check( nf90_get_var(ncid_in, varid, field_data, start = start), 'Read field' )
   else
     call check( nf90_get_var(ncid_in, varid, field_data), 'Read Field' )
   end if ! ndim > N

   end subroutine io_readfield_1d_real
   !******************************************************************
   ! 2D
   subroutine io_readfield_2d_real(field_name_in,time_dim,dims,field_data) bind(c)

   type(c_ptr), intent(in)      :: field_name_in
   integer(kind=c_int), value, intent(in) :: time_dim
   integer(kind=c_int), intent(in) :: dims(2)
   real(kind=c_real),intent(out) :: field_data(dims(1),dims(2))

   character(len=256)  :: field_name
   integer(kind=c_int) :: ndims
   integer(kind=c_int) :: varid
   integer(kind=c_int) :: this_dim

   call convert_c_string(field_name_in,field_name)
   this_dim = 2 

   call check( nf90_inq_varid(ncid_out, trim(field_name), varid), 'Get variable ID' )
   call check( nf90_inquire_variable(ncid_out,varid,ndims = ndims), 'Get variable DIMS' )
   if (ndims.gt.this_dim) then ! Then there is a time dimension, handle appropriately.
     start(:) = 1
     start(ndims) = time_dim
     call check( nf90_get_var(ncid_in, varid, field_data, start = start), 'Read field' )
   else
     call check( nf90_get_var(ncid_in, varid, field_data), 'Read Field' )
   end if ! ndim > N

   end subroutine io_readfield_2d_real
   !******************************************************************
   ! 3D
   subroutine io_readfield_3d_real(field_name_in,time_dim,dims,field_data) bind(c)

   type(c_ptr), intent(in)      :: field_name_in
   integer(kind=c_int), value, intent(in) :: time_dim
   integer(kind=c_int), intent(in) :: dims(3)
   real(kind=c_real),intent(out) :: field_data(dims(1),dims(2),dims(3))

   character(len=256)  :: field_name
   integer(kind=c_int) :: ndims
   integer(kind=c_int) :: varid
   integer(kind=c_int) :: this_dim

   call convert_c_string(field_name_in,field_name)
   this_dim = 3 

   call check( nf90_inq_varid(ncid_out, trim(field_name), varid), 'Get variable ID' )
   call check( nf90_inquire_variable(ncid_out,varid,ndims = ndims), 'Get variable DIMS' )
   if (ndims.gt.this_dim) then ! Then there is a time dimension, handle appropriately.
     start(:) = 1
     start(ndims) = time_dim
     call check( nf90_get_var(ncid_in, varid, field_data, start = start), 'Read field' )
   else
     call check( nf90_get_var(ncid_in, varid, field_data), 'Read Field' )
   end if ! ndim > N

   end subroutine io_readfield_3d_real
!======================================================================
   !******************************************************************
   ! 1D
   subroutine io_writefield_1d_real(field_name_in,time_dim,dims,field_data) bind(c)

   type(c_ptr), intent(in)      :: field_name_in
   integer(kind=c_int), value, intent(in) :: time_dim
   integer(kind=c_int), value, intent(in) :: dims
   real(kind=c_real),intent(in) :: field_data(dims)

   character(len=256)  :: field_name
   integer(kind=c_int) :: ndims
   integer(kind=c_int) :: varid
   integer(kind=c_int) :: this_dim

   call convert_c_string(field_name_in,field_name)
   this_dim = 1 

   call check( nf90_inq_varid(ncid_out, trim(field_name), varid), 'Get variable ID' )
   call check( nf90_inquire_variable(ncid_out,varid,ndims = ndims), 'Get variable DIMS' )
   if (ndims.gt.this_dim) then ! Then there is a time dimension, handle appropriately.
     start(:) = 1
     start(ndims) = time_dim
     call check( nf90_put_var(ncid_out, varid, field_data, start = start), 'Write field' )
   else
     call check( nf90_put_var(ncid_out, varid, field_data), 'Write Field' )
   end if ! ndim > N

   end subroutine io_writefield_1d_real
   !******************************************************************
   ! 2D
   subroutine io_writefield_2d_real(field_name_in,time_dim,dims,field_data) bind(c)

   type(c_ptr), intent(in)         :: field_name_in
   integer(kind=c_int), value, intent(in) :: time_dim
   integer(kind=c_int), intent(in) :: dims(2)
   real(kind=c_real),intent(in)    :: field_data(dims(1),dims(2))

   character(len=256)  :: field_name
   integer(kind=c_int) :: ndims
   integer(kind=c_int) :: varid
   integer(kind=c_int) :: this_dim

   call convert_c_string(field_name_in,field_name)
   this_dim = 2

   call check( nf90_inq_varid(ncid_out, trim(field_name), varid), 'Get variable ID' )
   call check( nf90_inquire_variable(ncid_out,varid,ndims = ndims), 'Get variable DIMS' )
   if (ndims.gt.this_dim) then ! Then there is a time dimension, handle appropriately.
     start(:) = 1
     start(ndims) = time_dim
     call check( nf90_put_var(ncid_out, varid, field_data, start = start), 'Write field' )
   else
     call check( nf90_put_var(ncid_out, varid, field_data), 'Write Field' )
   end if ! ndim > N

   end subroutine io_writefield_2d_real
   !******************************************************************
   ! 3D
   subroutine io_writefield_3d_real(field_name_in,time_dim,dims,field_data) bind(c)

   type(c_ptr), intent(in)         :: field_name_in
   integer(kind=c_int), value, intent(in) :: time_dim
   integer(kind=c_int), intent(in) :: dims(3)
   real(kind=c_real),intent(in)    :: field_data(dims(1),dims(2),dims(3))

   character(len=256)  :: field_name
   integer(kind=c_int) :: ndims
   integer(kind=c_int) :: varid
   integer(kind=c_int) :: this_dim

   call convert_c_string(field_name_in,field_name)
   this_dim = 3

   call check( nf90_inq_varid(ncid_out, trim(field_name), varid), 'Get variable ID' )
   call check( nf90_inquire_variable(ncid_out,varid,ndims = ndims), 'Get variable DIMS' )
   if (ndims.gt.this_dim) then ! Then there is a time dimension, handle appropriately.
     start(:) = 1
     start(ndims) = time_dim
     call check( nf90_put_var(ncid_out, varid, field_data, start = start), 'Write field' )
   else
     call check( nf90_put_var(ncid_out, varid, field_data), 'Write Field' )
   end if ! ndim > N

   end subroutine io_writefield_3d_real
!======================================================================
  subroutine convert_c_string(c_string_ptr,f_string)
  ! Purpose: To convert a c_string pointer to the proper fortran string format.
    type(c_ptr), intent(in) :: c_string_ptr
    character(len=256), intent(out) :: f_string
    character(len=256), pointer :: temp_string
    integer :: str_len

    call c_f_pointer(c_string_ptr,temp_string)
    str_len = index(temp_string, C_NULL_CHAR) - 1
    f_string = trim(temp_string(1:str_len))
    
    return
  end subroutine convert_c_string
!======================================================================
  subroutine convert_c_strarray(c_len,c_string_ptr,f_array)
  ! Purpose: To convert an array of c_string pointers to an array of fortran
  ! formatted strings.
    integer(kind=c_int), value, intent(in) :: c_len
    type(c_ptr), intent(in) :: c_string_ptr(c_len)
    character(len=256), dimension(c_len), intent(out) :: f_array
  
    integer :: i

    do i = 1,c_len
      call convert_c_string(c_string_ptr(i),f_array(i))
    end do

  return
  end subroutine convert_c_strarray
!======================================================================
  subroutine check(status, head)
  ! Purpose: Wrapper for standard netCDF commands with more detailed error
  ! messaging.
    integer(kind=c_int), intent ( in) :: status
    character(kind=c_char,len=*), intent(in) :: head

    character(kind=c_char,len=100) :: errmsg
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      errmsg = trim("Stopped at "//head)
      print *, errmsg
      stop 
    end if
  end subroutine check  
!======================================================================

end module simple_io_mod
