ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

include ../../make.inp

include ${OCCA_DIR}/scripts/Makefile

# define variables
HDRDIR = ../../
OGSDIR  = ../gatherScatter
GSDIR  = ../../3rdParty/gslib
HYPREDIR = ./hypre

# set options for this machine
# specify which compilers to use for c, fortran and linking
CC	= mpicc
CXX	= mpic++
LD	= mpic++

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I. -I./include/ $(compilerFlags) $(flags) -I$(HDRDIR)/include/ -I$(OGSDIR) -I$(GSDIR)/src -D DPARALMOND='"${CURDIR}"' $(LIBP_OPT_FLAGS)

# link flags to be used
LDFLAGS	= $(compilerFlags) $(flags) 

# libraries to be linked in
LIBS	=   -L$(OCCA_DIR)/lib  $(links) -L$(OGSDIR) -logs -L$(GSDIR)/lib -lgs -L$(HYPREDIR)/lib -lHYPRE \
	$(links) -L../../3rdParty/BlasLapack -lBlasLapack -lgfortran

INCLUDES = $(includes)
DEPS = $(INCLUDES) \
$(OGSDIR)/ogs.hpp \
$(HDRDIR)/include/types.h

# types of files we are going to construct rules for
.SUFFIXES: .c .cpp

# list of objects to be compiled
OBJS = \
./src/agmgLevel.o \
./src/agmgSmoother.o \
./src/coarseSolver.o \
./src/kernels.o \
./src/level.o \
./src/matrix.o \
./src/multigrid.o \
./src/parAlmond.o \
./src/pcg.o \
./src/pgmres.o \
./src/solver.o \
./src/SpMV.o \
./src/utils.o \
./src/vector.o \
./src/agmgSetup/agmgSetup.o \
./src/agmgSetup/constructProlongation.o \
./src/agmgSetup/formAggregates.o \
./src/agmgSetup/galerkinProd.o \
./src/agmgSetup/strongGraph.o \
./src/agmgSetup/transpose.o \
$(HYPREDIR)/hypre.o

ifeq ($(OS),Windows_NT)
    detected_OS := Windows
else
    detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
endif

SONAME=
WHOLE=
NOWHOLE=

ifeq ($(detected_OS),Linux)
	SONAME=-soname
	WHOLE=--whole-archive
	NOWHOLE=--no-whole-archive
	SHARED=-shared
	EXT=so
endif
ifeq ($(detected_OS),Darwin)
	SONAME=-install_name
	WHOLE=-all_load
	NOWHOLE=-noall_load
	SHARED=-dynamiclib -undefined dynamic_lookup
	EXT=dylib
endif

# rule for .c files
.c.o: $(DEPS) 
	$(CC) $(CFLAGS) -I$(HYPREDIR) -I$(HYPREDIR)/include -o $*.o -c $*.c $(paths)

.cpp.o: $(DEPS) 
	$(CXX) $(CFLAGS) -I$(HYPREDIR) -I$(HYPREDIR)/include -o $*.o -c $*.cpp $(paths)

.PHONY: all clean hypre

all: lib

hypre:
	cd hypre && env CC="$(CC)" NEKRS_SRC_THIRD_PARTY_DIR="${NEKRS_SRC_THIRD_PARTY_DIR}" \
	  ./install

lib: hypre
	@$(MAKE) --no-print-directory libParAlmond

libParAlmond: hypre  $(OBJS)
	$(CXX) $(SHARED) -o libparAlmond.$(EXT) $(OBJS) -L$(HYPREDIR)/lib -lHYPRE \
	    -Wl,$(SONAME),$(NEKRS_INSTALL_DIR)/parAlmond/libparAlmond.$(EXT)

clean:
	rm -f libparAlmond.a
	rm -f ./src/*.o
	rm -f ./src/agmgSetup/*.o
