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
  subroutine simple_io_get_field_size(fieldname_in,ncid,ndims,ldims,ndim_local) bind(c)

    type(c_ptr), intent(in)          :: fieldname_in
    integer(kind=c_int), value, intent(in)  :: ncid
    integer(kind=c_int), value, intent(in)  :: ndims
    integer(kind=c_int), intent(out) :: ldims(ndims)
    integer(kind=c_int), intent(out) :: ndim_local

    character(len=256) :: fieldname, dimname
    integer :: varid, dimids(ndims),ncid_loc
    integer :: n
  
    call convert_c_string(fieldname_in,fieldname)
    if (ncid.eq.0) then
      ncid_loc = ncid_in
    else
      ncid_loc = ncid_out
    end if
    
    call check( nf90_inq_varid(ncid_loc, fieldname, varid), "get VAR Id" )
    call check( nf90_inquire_variable(ncid_loc, varid, dimids=dimids, ndims=ndim_local), "get VAR-Dim Ids" )

    ldims(:) = 0
    do n = 1,ndim_local
      call check( nf90_inquire_dimension(ncid_loc, dimids(n), len=ldims(n), name=dimname), "get Var-Dim Len's" )
      if (trim(dimname)=='Time'.or.trim(dimname)=='time') then
        ldims(n) = -999
      end if
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
   ! Purpose: To read an input field from the input netCDF file and put it into
   ! the proper C++ format.
   ! Note that in order to grab a field from netCDF we need a local fortran
   ! array of at least the right dimension.  Since we can not declare a local
   ! array with arbitrary dimensionality, we assume a max dimension size of 6
   ! here.  At which point we only allocate for the total number of dimensions
   ! in the field.  Given that this is a simple interface that will be replaced
   ! it didn't make sense to spend too much time coming up with a more
   ! sophisticated way to handle this issue.
 
   type(c_ptr), intent(in)                :: field_name_in  ! Field name in C++ format
   integer(kind=c_int), intent(in)        :: time_dim(2)    ! Time dimension information which will be used for fields with a time dimension
   integer(kind=c_int), value, intent(in) :: flen           ! Full length of field in 1-D, i.e. flen=prod(field_dim_rng)
   real(kind=c_real),intent(out)          :: field_data(flen) ! 1-D vector to store field information

   real(kind=c_real), allocatable :: local_data(:,:,:,:,:,:)
   character(len=256)  :: field_name   ! Field name in fortran format
   integer             :: ndims, ndims_spatial,  varid ! Variable id and number of dimensions for variable
   integer             :: dimids(6), dimlens(6), start(6) ! Storage of dimension ids, lengths and starting point for read.  Note that start is needed to deal with variables that have a time dimension
   integer(kind=c_int) :: dimprod ! Check to make sure that field_data is long enough to fit all data.
   character(len=256)  :: dimname
   integer(kind=c_int) :: ii

   ! Convert field name to fortran string
   call convert_c_string(field_name_in,field_name)

   ! Retrieve the variable and dimension ids for this field
   call check( nf90_inq_varid(ncid_in, trim(field_name), varid), 'Get variable ID' )
   call check( nf90_inquire_variable(ncid_in,varid,ndims = ndims,dimids=dimids), 'Get variable DIMS' )
   
   ! Gather the max range for each dimension.
   dimprod = 1 ! Check that flen and dim lens match
   dimlens(:) = 0
   ndims_spatial = ndims
   do ii = 1,ndims
     call check( nf90_inquire_dimension(ncid_in, dimids(ii), len=dimlens(ii),name=dimname),  'Get dim lens' ) 
      if (trim(dimname).ne.'Time'.and.trim(dimname).ne.'time') then
        dimprod = dimprod*dimlens(ii)
       else
        ndims_spatial = ndims_spatial-1
      end if
   end do
   ndims = ndims_spatial
   if (dimprod.ne.flen) then
     print *, 'Field len = ', flen, ' and total dimension len = ',dimprod, ' does not match, exiting...'
     stop
   end if

   ! Allocate local data array to retrieve variable data
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

   ! Determine the starting position for all dimensions (important for time dimension)
   start(:) = 1
   start(time_dim(1)) = time_dim(2)
   ! Retrieve variable data from file and store in local data array.
   call check( nf90_get_var(ncid_in,varid, local_data, start=start), 'Get variable info')
   ! Reshape local_data array into a 1D vector to be passed out to C++
   field_data = reshape(local_data, (/ flen /))

   deallocate(local_data)

   return
   end subroutine io_readfield_real
!======================================================================
   subroutine io_writefield_real(field_name_in, time_dim, flen, field_data) bind(c)
   ! Purpose: To write an output field to the output netCDF file and put it into
   ! the proper C++ format.
   ! Note that in order to grab a field from netCDF we need a local fortran
   ! array of at least the right dimension.  Since we can not declare a local
   ! array with arbitrary dimensionality, we assume a max dimension size of 6
   ! here.  At which point we only allocate for the total number of dimensions
   ! in the field.  Given that this is a simple interface that will be replaced
   ! it didn't make sense to spend too much time coming up with a more
   ! sophisticated way to handle this issue.
 
   type(c_ptr), intent(in)                :: field_name_in  ! Field name in C++ format
   integer(kind=c_int), intent(in)        :: time_dim(2)    ! Time dimension information which will be used for fields with a time dimension
   integer(kind=c_int), value, intent(in) :: flen           ! Full length of field in 1-D, i.e. flen=prod(field_dim_rng)
   real(kind=c_real),intent(in)           :: field_data(flen) ! 1-D vector to store field information

   real(kind=c_real), allocatable :: local_data(:,:,:,:,:,:)
   character(len=256)  :: field_name   ! Field name in fortran format
   integer             :: ndims, ndims_spatial, varid ! Variable id and number of dimensions for variable
   integer             :: dimids(6), dimlens(6), start(6) ! Storage of dimension ids, lengths and starting point for read.  Note that start is needed to deal with variables that have a time dimension
   integer(kind=c_int) :: dimprod ! Check to make sure that field_data is long enough to fit all data.
   character(len=256)  :: dimname
   integer(kind=c_int) :: ii

   ! Convert field name to fortran string
   call convert_c_string(field_name_in,field_name)

   ! Retrieve the variable and dimension ids for this field
   call check( nf90_inq_varid(ncid_out, trim(field_name), varid), 'Get variable ID' )
   call check( nf90_inquire_variable(ncid_out,varid,ndims = ndims,dimids=dimids), 'Get variable DIMS' )
   
   ! Gather the max range for each dimension.
   dimprod = 1 ! Check that flen and dim lens match
   dimlens(:) = 0
   ndims_spatial = ndims
   do ii = 1,ndims
     call check( nf90_inquire_dimension(ncid_out, dimids(ii), len=dimlens(ii), name=dimname),  'Get dim lens' ) 
     if (trim(dimname).ne.'Time'.and.trim(dimname).ne.'time') then
       dimprod = dimprod*dimlens(ii)
       else
        ndims_spatial = ndims_spatial-1
     end if
   end do
   ndims = ndims_spatial
   if (dimprod.ne.flen) then
     print *, 'Field len = ', flen, ' and total dimension len = ',dimprod, ' does not match, exiting...'
     stop
   end if

   ! Allocate local data array to retrieve variable data
   select case (ndims)
     case (1)
       allocate(local_data(dimlens(1),1,1,1,1,1))
       local_data(:,1,1,1,1,1) = field_data
     case (2)
       allocate(local_data(dimlens(1),dimlens(2),1,1,1,1))
       local_data(:,:,1,1,1,1) = reshape(field_data, dimlens(:2))
     case (3)
       allocate(local_data(dimlens(1),dimlens(2),dimlens(3),1,1,1))
       local_data(:,:,:,1,1,1) = reshape(field_data, dimlens(:3))
     case (4)
       allocate(local_data(dimlens(1),dimlens(2),dimlens(3),dimlens(4),1,1))
       local_data(:,:,:,:,1,1) = reshape(field_data, dimlens(:4))
     case (5)
       allocate(local_data(dimlens(1),dimlens(2),dimlens(3),dimlens(4),dimlens(5),1))
       local_data(:,:,:,:,:,1) = reshape(field_data, dimlens(:5))
     case (6)
       allocate(local_data(dimlens(1),dimlens(2),dimlens(3),dimlens(4),dimlens(5),dimlens(6)))
       local_data(:,:,:,:,:,:) = reshape(field_data, dimlens(:6))
     case default
       print *, "WRITEFIELD: ndims > 6 not supported"
       stop
   end select

   ! Determine the starting position for all dimensions (important for time dimension)
   start(:) = 1
   start(time_dim(1)) = time_dim(2)
   ! Retrieve variable data from file and store in local data array.
   call check( nf90_put_var(ncid_out, varid, local_data, start = start), 'Write field' )

   deallocate(local_data)

   return
   end subroutine io_writefield_real
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
