#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module simple_io_mod
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
    ! Simple routine to open a netcdf file for reading input.
    type(c_ptr), intent(in) :: filename_in
    
    character(len=256) :: filename

    call convert_c_string(filename_in,filename)
    call check( nf90_open(filename, nf90_nowrite, ncid_in), "Open input file" )
    isopen_in = .true.

    return
  end subroutine simple_in_init
!======================================================================
  subroutine simple_in_finalize() bind(c)
    ! Simple routine to close a netcdf file from which input has been read.

    call check( nf90_close(ncid_in), "Close input file" )
    isopen_in = .false.

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
    ! Simple routine to open (and replace if needed) a netcdf output file.
    ! Register all the potential dimensions
  
    type(c_ptr), intent(in) :: filename_in
    integer(kind=c_int), value, intent(in)  :: ndims  ! Number of potential dimensions
    type(c_ptr), intent(in) :: dimnames_in(ndims)  ! Name for each dimension
    integer(kind=c_int), target, intent(in)  :: dimrng(ndims) ! Max range for each dimension
  
    character(len=256) :: dimnames(ndims)
    character(len=256) :: filename
  
    integer(kind=c_int) :: dimids(ndims) ! The id for each dimension
    
    integer(kind=c_int) :: i
  
    call convert_c_string(filename_in,filename)
    call convert_c_strarray(ndims,dimnames_in,dimnames)
  
    call check( nf90_create(trim(filename), nf90_clobber, ncid_out), 'init1' )
    isopen_out = .true.
  
    do i = 1,ndims
      call check( nf90_def_dim(ncid_out, trim(dimnames(i)), dimrng(i), dimids(i)),'add dimensions' )
    end do
  
    return
  end subroutine simple_out_init1
!======================================================================
  subroutine simple_out_init2() bind(c)
  
    allocate(start(max_ndims))
  
    call check( nf90_enddef(ncid_out), 'init2' )
    isopen_out = .false.
  
  end subroutine simple_out_init2
!======================================================================
  subroutine simple_out_finalize() bind(c)
  
    call check( nf90_close(ncid_out), 'finalize' )
  
  end subroutine simple_out_finalize
!======================================================================
  subroutine simple_out_regfield(field_name_in,field_type,ndim,field_dim_in,units_in) bind(c)
  
    type(c_ptr),         intent(in)           :: field_name_in
    integer(kind=c_int), value, intent(in)    :: field_type
    integer(kind=c_int), value, intent(in)    :: ndim
    type(c_ptr),         intent(in)           :: field_dim_in(ndim)
    type(c_ptr),         intent(in), optional :: units_in
  
    character(len=256) :: field_name
    character(len=256) :: field_dim(ndim)
    character(len=256) :: units
  
    integer(kind=c_int) :: varid
    integer(kind=c_int) :: dimids(ndim)
    character(kind=c_char,len=100) :: check_case
    integer(kind=c_int) :: i
  
    if (.not.isopen_out) then
      print *, 'netCDF definition is closed, cannot register new field', field_name
      stop
    end if
  
    call convert_c_string(field_name_in,field_name)
    call convert_c_strarray(ndim,field_dim_in,field_dim)
    if (present(units_in)) call convert_c_string(units_in,units)
  
    do i = 1,ndim
      call check( nf90_inq_dimid(ncid_out,trim(field_dim(i)),dimids(i)), 'Reg: get dim ids' )
    end do
  
    check_case = 'Add variable '//trim(field_name)
    call check( nf90_def_var(ncid_out,trim(field_name), field_type, dimids, varid), check_case )
    check_case = 'Add variable '//trim(field_name)//', units'
    if (present(units_in)) then
      call check( nf90_put_att(ncid_out, varid, "units", trim(units)), check_case )
    else
      call check( nf90_put_att(ncid_out, varid, "units", "unspecified"), check_case )
    end if 
  
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
