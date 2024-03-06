# Makefile -- Use this to build on *NIX systems.

# Options set on command line.
debug      = 1
mpi        = not-set
shared     = not-set
precision  = not-set
verbose    = not-set
prefix     = not-set
sanitize   = not-set
CC         = not-set
CXX        = not-set
FC         = not-set
travis     = not-set
# This proxies everything to the builddir cmake.
netcdfsys = not-set

cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)

BUILDDIR := build/$(systype)-$(cputype)
CONFIG_FLAGS = -DUNIX=1 -Wno-dev

# Process configuration options.

# Travis-CI build
ifeq ($(travis), not-set)
  CONFIG_FLAGS += -DTRAVIS_CI=0
else
  CONFIG_FLAGS += -DTRAVIS_CI=1
endif

# Verbose builds?
ifeq ($(verbose), 1)
  CONFIG_FLAGS += -DCMAKE_VERBOSE_MAKEFILE=1
endif

# MPI
ifeq ($(mpi), 1)
  BUILDDIR := ${BUILDDIR}-mpi
  CC = mpicc
  CXX = mpicxx
  FC = mpif90
  CONFIG_FLAGS += -DHAVE_MPI=1
else
  ifeq ($(CC), not-set)
    CC  = icc
  endif
  ifeq ($(CXX), not-set)
    CXX = icpc
  endif
  ifeq ($(FC), not-set)
    FC = ifort
  endif
  CONFIG_FLAGS += -DHAVE_MPI=0
endif

# Shared libs?
ifeq ($(shared), 1)
  BUILDDIR := ${BUILDDIR}-shared
  CONFIG_FLAGS += -DBUILD_SHARED_LIBS=ON
else
  BUILDDIR := ${BUILDDIR}-static
  CONFIG_FLAGS += -DBUILD_SHARED_LIBS=OFF
endif

# Precision.
ifneq ($(precision), not-set)
  BUILDDIR := ${BUILDDIR}-double
  CONFIG_FLAGS += -DECOSYS_PRECISION=double
else
  BUILDDIR := ${BUILDDIR}-$(precision)
  CONFIG_FLAGS += -DECOSYS_PRECISION=$(precision)
endif

BUILDDIR := ${BUILDDIR}-`basename ${CC}`
CONFIG_FLAGS += -DCC=${CC} -DCXX=${CXX}
ifneq ($(FC), )
  CONFIG_FLAGS += -DFC=${FC}
endif

# Debugging symbols
ifeq ($(debug), not-set)
  BUILDDIR := ${BUILDDIR}-Release
  CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Release
else
  ifeq ($(debug), 0)
    BUILDDIR := ${BUILDDIR}-Release
    CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Release
  else
    BUILDDIR := ${BUILDDIR}-Debug
    CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Debug
  endif
endif

# Installation prefix.
ifeq ($(prefix), not-set)
  prefix = $(CURDIR)/local
endif
CONFIG_FLAGS += -DCMAKE_INSTALL_PREFIX:PATH=$(prefix)

# Special considerations for specific systems.
ifeq ($(systype), Darwin)
  CONFIG_FLAGS += -DAPPLE=1
else
  ifeq ($(systype), Linux)
    CONFIG_FLAGS += -DLINUX=1
  endif
endif

# Address sanitizer.
ifeq ($(sanitize), 1)
  BUILDDIR := ${BUILDDIR}-AddressSanitizer
  CONFIG_FLAGS += -DADDRESS_SANITIZER=1
endif

# netcdf

#ifeq ($(netcdfsys), not-set)
#  NETCDF_FFLAGS =""
#  NETCDF_FLIBS =""
#else
#  NETCDF_FFLAGS = $(shell nc-config --prefix)/include
#  NETCDF_FLIBS = $(shell nc-config --flibs)
#endif

#CONFIG_FLAGS += -DTPL_NETCDF_INCLUDE_DIRS=" $(NETCDF_FFLAGS)"
#CONFIG_FLAGS += -DTPL_NETCDF_LIBRARIES="$(NETCDF_FLIBS)"

define run-config
@mkdir -p $(BUILDDIR)
@cd $(BUILDDIR) && cmake $(CURDIR) $(CONFIG_FLAGS)
endef

all:
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more INSTALL; \
	else \
		$(MAKE) -C $(BUILDDIR) $@ --no-print-directory $(MAKEFLAGS); \
	fi

install: all
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more INSTALL; \
	else \
		$(MAKE) -C $(BUILDDIR) $@ --no-print-directory $(MAKEFLAGS); \
	fi

test: install
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more INSTALL; \
	else \
		$(MAKE) -C $(BUILDDIR) $@ --no-print-directory $(MAKEFLAGS); \
		$(MAKE) -C regression-tests $@ --no-print-directory $(MAKEFLAGS); \
	fi

clean:
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more INSTALL; \
	else \
		$(MAKE) -C $(BUILDDIR) $@ --no-print-directory $(MAKEFLAGS); \
	fi

config: distclean
	$(run-config)
	@export ENVCC=$(CC); export ENVCXX=$(CXX); export ENVFC=$(FC)

distclean:
	@rm -rf $(BUILDDIR)
	@rm -rf ./local

stats:
	@python tools/gather_stats.py

prepend-license:
	@python tools/prepend_license.py

ctags-emacs :
	@ctags -e -f ETAGS -R --exclude=.git --exclude=build

#dist:
#	utils/mkdist.sh $(PKGNAME)

.PHONY: config distclean all clean install uninstall test
