#
# BUG?!?! 
#
# 1) if you get random arlib/link failures rebuild MARBL source in its dir. This can 
# happen after OSX or Matlab verssion updates and so on.
#
# 
# This Makefile uses gfortran to build the MARBL library with the -fPIC flag
# The mex target then uses these object files when building the marbl interface
# for use with mex / Matlab
#
# Assumes that MARBL is available in the directory marbl parallel to
# marbl-interface/. If this is not the case, run
# $ make MARBL_ROOT=[marbl directory]
# 
# NOTE: some magic in Matlab MEX changes MEX functions to have these 64-bit,
# platform-specific file extensions:
# 
# Linux® — .mexa64
# Apple macOS — .mexmaci64
# Microsoft® Windows® — .mexw64

# For my Mac
# 		MARBL_ROOT=../../MARBL-stable_2020_07_27
# This is for Green Planet
# 		MARBL_ROOT=/nfspool-0/home/jjbecker/Desktop/MARBL/MARBL-stable_2020_07_27

MARBL_ROOT=../../MARBL-stable_2020_07_27
ifneq ($(shell uname -s),Darwin)
    MARBL_ROOT=../../../../Desktop/MARBL/MARBL-stable_2020_07_27
endif
$(info $$MARBL_ROOT is [${MARBL_ROOT}])


MARBL_ROOT_FULL=$(realpath $(MARBL_ROOT))
LIB_DIR=$(realpath .)/marbl_lib
MARBL_LIB=marbl_lib/libmarbl.a
SO_INTERFACE=marbl_lib/marbl_interface.so
INTERFACE_SRC=marbl_interface_wrapper_mod.F90
INTERFACE_OBJ=$(INTERFACE_SRC:.F90=.o)
MEX_INTERFACE=$(INTERFACE_SRC:.F90=.mexa64)
MEX_DRIVER_SRC=mex_marbl_driver.F90
MEX_DRIVER=$(MEX_DRIVER_SRC:.F90=.mexa64)

INC_DIR = marbl_include
dummy_build_marbl_include := $(shell mkdir -p $(INC_DIR))

dummy_build_marbl_lib := $(shell mkdir -p $(LIB_DIR))


# This is for gfortran 8.3 
FFLAGS += -fPIC -Wall -fmax-stack-var-size=8388607
# This is for gfortran 10.2 
# FFLAGS += $(FFLAGS) -fallow-invalid-boz

myLDFLAGS=""
# myLDFLAGS="-Wl,-no_compact_unwind -fallow-invalid-boz -Wall -fmax-stack-var-size=8388607"


#############
#  TARGETS  #
#############

# By default, build the mex driver included as an example
all: mex

# The mex driver should be rebuilt if the mex interface changes or if the
# driver source code changes
#
# FIXME: huge strings for MAX_TMP_CNT and MAX_LOG_ENTRY_CNT cause crazy 
# big variables on the stack which cause bugs if multithreaded, hence 
# replace default make variable size with 8 MB-1, 
#   -q because default is 64K -1
#
$(MEX_DRIVER): $(MEX_INTERFACE) $(MEX_DRIVER_SRC)
	$(info $$MARBL_ROOT_FULL is [${MARBL_ROOT_FULL}])
	mex FOPTIMFLAGS='-O3 $(FFLAGS)' -Imarbl_include $(MEX_DRIVER_SRC) $(INTERFACE_SRC) marbl_include/*.o 

# Shortcut for building the mex driver: "$ make mex"
.PHONY: mex
mex: $(MEX_DRIVER)

# The mex interface should be rebuilt if the MARBL library changes or if
# the interface source code changes
$(MEX_INTERFACE): $(MARBL_LIB) $(INTERFACE_SRC)
	$(info $$MARBL_ROOT_FULL is [${MARBL_ROOT_FULL}])
	mex FOPTIMFLAGS='-O3 $(FFLAGS)' -Imarbl_include $(INTERFACE_SRC) -Lmarbl_lib -lmarbl -c 

$(SO_INTERFACE): $(MARBL_LIB) $(INTERFACE_SRC)
	$(info $$MARBL_ROOT_FULL is [${MARBL_ROOT_FULL}])
	gfortran -Imarbl_include $(FFLAGS) -c $(INTERFACE_SRC) 
	gfortran -Lmarbl_lib -shared -o $(SO_INTERFACE) $(INTERFACE_OBJ) -lmarbl

# Short target for building the .so target
libso: $(SO_INTERFACE)

# The marbl library should be rebuilt if any of the MARBL fortran files change
$(MARBL_LIB): $(wildcard $(MARBL_ROOT_FULL)/src/*.F90)
	$(info $$MARBL_ROOT_FULL is [${MARBL_ROOT_FULL}])
	cd marbl_include ; pwd ; ls ../marbl_lib ; \
        $(MAKE) -f $(MARBL_ROOT_FULL)/src/Makefile \
            FC=gfortran FCFLAGS="-O3 $(FFLAGS)" \
                LDFLAGS=$(myLDFLAGS) \
            USE_DEPS=TRUE \
            OBJ_DIR=. INC_DIR=. LIB_DIR=../marbl_lib ../$(MARBL_LIB) ; \
    cd ..

# Here's an easy way to build just the MARBL library: "$ make lib"
# (Instead of "$ make marbl_lib/libmarbl.a")
.PHONY: lib
lib: $(MARBL_LIB)

# Clean up just the interface files with "$ make clean"
.PHONY: clean
clean:
	$(info $$MARBL_ROOT_FULL is [${MARBL_ROOT_FULL}])
	rm -f *.mod *.mexmaci64 *.o *.pyc $(SO_INTERFACE)

# Or use "$ make allclean" to clean up the interface files and the MARBL library
.PHONY: allclean
allclean: clean
	$(info $$MARBL_ROOT_FULL is [${MARBL_ROOT_FULL}])
	rm -f marbl_include/*.mod marbl_include/*.o marbl_include/*.d $(MARBL_LIB)
	rmdir marbl_include marbl_lib

