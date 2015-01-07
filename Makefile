
# RDIAG Toolkit Makefile principal

ifeq "$(BASE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

SHELL   = /bin/bash

# Location of the development versions of the src, bin/$(BASE_ARCH),
# lib/$(EC_ARCH) and man/pdoc directory trees

DIAGNOSTIQUE = $(CURDIR)
MAKE    = make DIAGNOSTIQUE=$(DIAGNOSTIQUE)

# Destination that the current working binaries and libraries
# will be exported to, as well as where to find the EXTRAS

DESTINATION = $(DIAGNOSTIQUE)/..

# NetCDF v_3.6 and Udunits v_1.2 libraries (the EXTRAS)

EXTRAS  = $(DESTINATION)/extras
EXTLIB  = $(EXTRAS)/NetcdfUdunits/$(EC_ARCH)/lib

# Directories used/created by this Makefile

BINDIR  = $(DIAGNOSTIQUE)/bin/$(BASE_ARCH)
INCLUDE = $(DIAGNOSTIQUE)/include
LIBDIR  = $(DIAGNOSTIQUE)/lib/$(EC_ARCH)
MANDIR  = $(DIAGNOSTIQUE)/man/pdoc
MODDIR  = $(INCLUDE)/$(EC_ARCH)
SUBDIR  = $(DIAGNOSTIQUE)/src/lssub/sources$(STD)/$(EC_ARCH)

# Include (very) old NCAR graphics (default=no) ?

NOPLOT  = -DNOPLOT
GRAFLIB =

# Binaires pre-compiles

FIXES   =

# (obsolete) WEB Host Server for the documentation

#HOSTWEB = notos
#DIAGWEB = /data/armnraid1/www/mrb/externe/si/eng/si/utilities/r.diag

# RMN and Vgrid_Descriptor library names

RMNLIB  = rmn_015
VGDLIB  = descrip

# DDFUN90, NetCDF3 and UdUnits1 library names

DDFUN90  = ddfun90
lNetCDF  = netcdff netcdf
UDUNITS  = udunits

DIAG_VERSION = 6.3.0
CONV_VERSION = 2.2.0

STD     = 98

default: allbin

allbin: initial rdiag cdf2conv 

all: allbin document

export:
	/bin/mkdir -p $(DESTINATION)/lib $(DESTINATION)/bin ;\
	rsync -av $(DIAGNOSTIQUE)/lib/$(EC_ARCH) $(DESTINATION)/lib/$(BASE_ARCH) ;\
	rsync -av $(DIAGNOSTIQUE)/bin/$(BASE_ARCH) $(DESTINATION)/bin ;\
	rsync -av $(DIAGNOSTIQUE)/man $(DESTINATION)

# Ensure initial setup is done

initial:
	/bin/mkdir -p $(BINDIR) $(LIBDIR) $(MANDIR) $(MODDIR) $(SUBDIR)
	s.locate --lib $(VGDLIB) 1> /dev/null || { echo -e PLS execute \". s.ssmuse.dot vgriddesc\" \n ; false ; }
	if [[ ! -f $(EXTLIB)/libnetcdf.a ]]; then cd $(EXTRAS) ; make all ; fi
	if [[ ! -f $(LIBDIR)/libddfun90.a || -z "$(DDFUN90)" ]]; then \
	cd $(DIAGNOSTIQUE)/src/extras/ddfun90 ; $(MAKE) ; fi
	if [[ ! -x $(BINDIR)/r.echo ]]; then cd $(DIAGNOSTIQUE)/src/extras/tools ; $(MAKE) ; fi
	if [[ ! -f $(LIBDIR)/program_version.o ]]; then cd $(LIBDIR) ;\
	s.f77 -g -c ../../program_version.f ; fi

# RDIAG Diagnostic toolkit recipe

rdiag:
	echo "*** Making libdiag_sq98.a and libdiag_sq98_g.a ***" ;\
	cd $(DIAGNOSTIQUE)/src/lssub ; $(MAKE) VGDLIB=$(VGDLIB)
	echo "Making libprog_sq98.a" ;\
	cd $(DIAGNOSTIQUE)/src/lspgm ; $(MAKE)
	echo "*** Making executable r.diag ***" ;\
	cd $(DIAGNOSTIQUE)/src/lspgm ; $(MAKE) $(BASE_ARCH) OBJ=$(FIXES) \
	NOPLOT=$(NOPLOT) GRAFLIB=$(GRAFLIB) DDFUN90=$(DDFUN90) VGDLIB=$(VGDLIB) \
	DIAG_VERSION=$(DIAG_VERSION) RMNLIB=$(RMNLIB)

# NetCDF to/from ( CCCma or CMC/RPN) file format converter recipe

cdf2conv:
	echo "*** Making libcdf2ccc.a ***" ;\
	cd $(DIAGNOSTIQUE)/src/cdf2ccc ; $(MAKE) \
	EXTRAS=$(EXTRAS)/NetcdfUdunits/$(EC_ARCH)
	echo "*** Making executable cdf2ccc ***" ;\
	cd $(DIAGNOSTIQUE)/src/cdf2ccc ;\
	$(MAKE) cdf2rpn CONV_VERSION=$(CONV_VERSION) \
	RMNLIB=$(RMNLIB) VGDLIB=$(VGDLIB) OBJ=$(FIXES) \
	EXTRAS=$(EXTRAS)/NetcdfUdunits/$(EC_ARCH) DDFUN90=$(DDFUN90) \
	lNetCDF="$(lNetCDF)" UDUNITS=$(UDUNITS)

# Online documentation (which was originaly found in $ARMNLIB/man/pdoc) recipe

document:
	cd $(DIAGNOSTIQUE)/src/lspgm   ; $(MAKE) $@ ; $(MAKE) info.lspgm ;\
	cd $(DIAGNOSTIQUE)/src/lssub   ; $(MAKE) $@ ;\
	cd $(DIAGNOSTIQUE)/src/cdf2ccc ; $(MAKE) $@
