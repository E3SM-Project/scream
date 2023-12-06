# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file
. /usr/share/lmod/lmod/init/sh
module load python/3.9.12 git mkl/2022.1.0 intel-classic/2021.6.0-magic mvapich2/2.3.7 cmake/3.19.2 netcdf-fortran-parallel/4.6.0 netcdf-c-parallel/4.9.0 parallel-netcdf/1.12.3
export NETCDFROOT=/usr/tce/packages/netcdf-fortran/netcdf-fortran-4.6.0-mvapich2-2.3.7-intel-classic-2021.6.0/
export NETCDF_PATH=/usr/tce/packages/netcdf-fortran/netcdf-fortran-4.6.0-mvapich2-2.3.7-intel-classic-2021.6.0/
export PNETCDFROOT=/usr/tce/packages/parallel-netcdf/parallel-netcdf-1.12.3-mvapich2-2.3.7-intel-classic-2021.6.0/
export COMPILER=intel
export MPILIB=mpi-serial
export DEBUG=FALSE
export OS=LINUX
