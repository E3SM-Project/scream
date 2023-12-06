CFLAGS :=  -fp-model precise -std=gnu99 -O2 -debug minimal
CPPDEFS := $(CPPDEFS)  -DFORTRANUNDERSCORE -DNO_R16 -DCPRINTEL -DNO_SHR_VMATH -DCNL
CXXFLAGS :=  -fp-model source -O2
CXX_LDFLAGS :=  -cxxlib
CXX_LINKER := FORTRAN
FC_AUTO_R8 :=  -r8
FFLAGS :=  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source -O2 -debug minimal
FFLAGS_NOOPT :=  -O0
FIXEDFLAGS :=  -fixed -132
FREEFLAGS :=  -free
HAS_F2008_CONTIGUOUS := TRUE
KOKKOS_OPTIONS := --with-serial --ldflags='-L/usr/tce/packages/gcc/gcc-10.3.1-magic/lib/gcc/x86_64-redhat-linux/10/'
LDFLAGS :=  -L/usr/tce/packages/gcc/gcc-10.3.1-magic/lib/gcc/x86_64-redhat-linux/10/
MACRO_FILE := 
MPICC := mpicc
MPICXX := mpicxx
MPIFC := mpif90
MPI_LIB_NAME := mpich
MPI_PATH := /usr/tce/packages/mvapich2/mvapich2-2.3.7-intel-classic-2021.6.0/
NETCDF_PATH := /usr/tce/packages/netcdf-fortran/netcdf-fortran-4.6.0-mvapich2-2.3.7-intel-classic-2021.6.0/
PNETCDF_PATH := /usr/tce/packages/parallel-netcdf/parallel-netcdf-1.12.3-mvapich2-2.3.7-intel-classic-2021.6.0/
SCC := icc
SCXX := icpc
SFC := ifort
SLIBS := $(SLIBS)  -llapack -lblas -L/usr/tce/backend/installations/linux-rhel8-x86_64/intel-2021.6.0/netcdf-fortran-4.6.0-gbg3y7c2xjnda55wkwwzupgzgcmevm3w/lib -lnetcdff -lnetcdf -lnetcdf -lm
SUPPORTS_CXX := TRUE

