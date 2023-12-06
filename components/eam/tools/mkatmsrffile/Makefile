#-------------------------------------
# Makefile for mkatmsrffile tool 
#-------------------------------------

#-------------------------------------
#Instructions:
# Load a compiler and the associated 
# netcdf library. Update variables "FC"
# and "NETCDF_PATH" to add the compiler
# and path to netcdf home directory
#-------------------------------------

ifeq ($(FC),$(null))
   FC := ifort
endif
FC=ifort

FFLAGS := $(shell nf-config --flibs)
LDFLAGS := $(shell nf-config --flibs)

NETCDF_PATH := $(NETCDF_HOME)

EXE := mkatmsrffile

#---------------------------
# USER INPUT ENDS HERE
#---------------------------

#include string
NETCDF_ROOT=/usr/tce/backend/installations/linux-rhel8-x86_64/intel-2021.6.0/netcdf-fortran-4.6.0-gbg3y7c2xjnda55wkwwzupgzgcmevm3w
#INC_NETCDF := $(shell nf-config --includedir)
#LIB_NETCDF := $(shell nc-config --libdir)
#NETCDF_ROOT=$(dirname $(dirname $INC_NETCDF))

INC = -I$(NETCDF_ROOT)/include
LIB = -L$(NETCDF_ROOT)/lib -lnetcdf -lnetcdff

#(zhang73)
FFLAGS := -I$(INC_NETCDF) 

FFLAGS += -O2


$(EXE): mkatmsrffile.F90 
	$(FC) $(FFLAGS) $(INC)  $< -o $@  $(LIB)

clean:
	/bin/rm $(EXE)

