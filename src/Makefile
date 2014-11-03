
ifeq "$(BASE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

SHELL   = /bin/bash

NOPLOT  = -DNOPLOT

# Location of the src, bin/$(BASE_ARCH), lib/$(EC_ARCH) and man/pdoc directory trees

DIAGNOSTIQUE = $(CURDIR)/..
MAKE    = make DIAGNOSTIQUE=$(DIAGNOSTIQUE)

EXTRAS  = $(DIAGNOSTIQUE)/extras
EXTLIB  = $(EXTRAS)/NetcdfUdunits/$(EC_ARCH)/lib

INCLUDE = $(DIAGNOSTIQUE)/include
LIBDIR  = $(DIAGNOSTIQUE)/lib/$(EC_ARCH)
BINDIR  = $(DIAGNOSTIQUE)/bin/$(BASE_ARCH)
MANDIR  = $(DIAGNOSTIQUE)/man/pdoc
MODDIR  = $(INCLUDE)/$(EC_ARCH)

# Serveur WEB hote pour la documentation

#HOSTWEB = notos
#DIAGWEB = /data/armnraid1/www/mrb/externe/si/eng/si/utilities/r.diag

# RMN and Vgrid_Descriptor library names

RMNLIB  = rmn_015
VGDLIB  = descrip

# NetCDF3 and UdUnits1 library names

lNetCDF  = netcdff netcdf
UDUNITS  = udunits

DIAG_VERSION = 6.3.0
CONV_VERSION = 2.2.0

STD     = 98

default: allbin

allbin: initial rdiag cdf2conv 

all: allbin document

# Ensure initial setup is done

initial:
	/bin/mkdir -p lssub/sources$(STD)/$(EC_ARCH) $(LIBDIR) $(MANDIR) $(MODDIR)
	s.locate --lib $(VGDLIB) 1> /dev/null || { echo -e PLS execute \". s.ssmuse.dot vgriddesc\" \n ; false ; }
	if [[ ! -f $(EXTLIB)/libnetcdf.a ]]; then cd $(EXTRAS) ; make all ; fi
	if [[ ! -f $(LIBDIR)/libddfun90.a ]]; then cd $(CURDIR)/extras/ddfun90 ; make ; fi
	if [[ ! -x $(BINDIR)/r.echo ]]; then cd $(CURDIR)/extras/tools ; make ; fi
	if [[ ! -f $(LIBDIR)/program_version.o ]]; then cd $(LIBDIR) ;\
	s.f90 -g -c ../../program_version.f ; fi

# RDIAG Diagnostic toolkit recipe

rdiag:
	echo "*** Making libdiag_sq98.a and libdiag_sq98_g.a ***" ;\
	cd $(DIAGNOSTIQUE)/src/lssub ; $(MAKE) VGDLIB=$(VGDLIB)
	echo "Making libprog_sq98.a" ;\
	cd $(DIAGNOSTIQUE)/src/lspgm ; $(MAKE)
	echo "*** Making executable r.diag ***" ;\
	cd $(DIAGNOSTIQUE)/src/lspgm ; $(MAKE) $(BASE_ARCH) NOPLOT=$(NOPLOT) \
	DIAG_VERSION=$(DIAG_VERSION) VGDLIB=$(VGDLIB) RMNLIB=$(RMNLIB)

# NetCDF to/from ( CCCma or CMC/RPN) file format converter recipe

cdf2conv:
	echo "*** Making libcdf2ccc.a ***" ;\
	cd $(DIAGNOSTIQUE)/src/cdf2ccc ; $(MAKE)
	echo "*** Making executable cdf2ccc ***" ;\
	cd $(DIAGNOSTIQUE)/src/cdf2ccc ;\
	$(MAKE) cdf2rpn CONV_VERSION=$(CONV_VERSION) \
	RMNLIB=$(RMNLIB) VGDLIB=$(VGDLIB) \
	EXTRAS=$(EXTRAS)/NetcdfUdunits/$(EC_ARCH) \
	lNetCDF="$(lNetCDF)" UDUNITS=$(UDUNITS) ; cd ..

# Online documentation (which was originaly found in $ARMNLIB/man/pdoc) recipe

document:
	cd $(DIAGNOSTIQUE)/src/lspgm   ; $(MAKE) $@ ; $(MAKE) info.lspgm ;\
	cd $(DIAGNOSTIQUE)/src/lssub   ; $(MAKE) $@ ;\
	cd $(DIAGNOSTIQUE)/src/cdf2ccc ; $(MAKE) $@
