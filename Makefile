
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

# NetCDF v_3.6 and Udunits v_1.2 libraries used to be found 
# in the EXTRAS directory. The package now can either 1) use
# the netcdff SSM package to locate all the necessary shared
# libraries including MFV's udunits2f FORTRAN wrapper as well
# as udunits2 itself or 2) use the static versions of these
# same libraries as provided by the netcdff-4.4 SSM package.
# As of May 2017, the second option is used.
#EXTRAS  = $(DESTINATION)/extras
#EXTLIB  = $(EXTRAS)/NetcdfUdunits/$(EC_ARCH)/lib

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

#FIXES   = $(PWD)/lib/$(EC_ARCH)/stubs.o

# (obsolete) WEB Host Server for the documentation

HOSTWEB = pascal
DIAGWEB = public_html

# RMN and Vgrid_Descriptor library names

RMNLIB  = rmn_016.3
VGDLIB  = descrip

# DDFUN90, NetCDF4 and UdUnits2 library names

DDFUN90  = ddfun90

# Shared load via the netcdff SSM package
#lNetCDF  = netcdff
#UDUNITS  = udunits2f udunits2

# Static (_s) load via symlinks in the netcdff-4.4 SSM package
lNetCDF  = netcdff_s netcdf_s hdf5_hl_s hdf5_s dl sz_s z
UDUNITS  = udunits2f_s udunits2_s expat


DIAG_VERSION = 6.4.1
CONV_VERSION = 2.3.1

ENTETE  = 32
STD     = 98

default: allbin

allbin: initial_base initial_cdf rdiag cdf2conv 

all: allbin document

export:
	/bin/mkdir -p $(DESTINATION)/lib $(DESTINATION)/bin ;\
	rsync -av $(DIAGNOSTIQUE)/lib/$(EC_ARCH) $(DESTINATION)/lib/$(BASE_ARCH) ;\
	rsync -avH $(DIAGNOSTIQUE)/bin/$(BASE_ARCH) $(DESTINATION)/bin ;\
	rsync -lptgoDv $(DIAGNOSTIQUE)/bin/* $(DESTINATION)/bin ;\
	rsync -av $(DIAGNOSTIQUE)/man $(DESTINATION)

# Ensure initial setup is done

initial_base:
	/bin/mkdir -p $(BINDIR) $(LIBDIR) $(MANDIR) $(MODDIR) $(SUBDIR)
	s.locate --lib $(VGDLIB) 1> /dev/null || { echo -e "\nPLS execute \". s.ssmuse.dot vgriddesc\"\n" ; false ; }
	s.locate --lib netcdff_s 1> /dev/null || { echo -e "\nPLS execute \". s.ssmuse.dot netcdff-4.4\"\n" ; false ; }
#	s.locate --lib $(lNetCDF) 1> /dev/null || { echo -e "\nPLS execute \". s.ssmuse.dot netcdff\"\n" ; false ; }
#	if [[ ! -f $(EXTLIB)/libnetcdf.a ]]; then cd $(EXTRAS) ; make all ; fi
	if [[ ! -f $(LIBDIR)/libddfun90.a || -z "$(DDFUN90)" ]]; then \
	cd $(DIAGNOSTIQUE)/src/extras/ddfun90 ; $(MAKE) RMNLIB=$(RMNLIB) ; fi
	if [[ ! -x $(BINDIR)/r.echo ]]; then cd $(DIAGNOSTIQUE)/src/extras/tools ; $(MAKE) ; fi
	if [[ ! -f $(LIBDIR)/program_version.o ]]; then cd $(LIBDIR) ;\
	s.f77 -g -c ../../program_version.f ; fi
#	if [[ ! -f $(LIBDIR)/crc32.o ]]; then cd $(LIBDIR) ;\
#	s.cc -g -c ../../crc32.c ; fi

initial_cdf:
#	if [[ ! -f $(EXTLIB)/libnetcdf.a ]]; then cd $(EXTRAS) ; make all ; fi
#	s.locate --lib $(lNetCDF) 1> /dev/null || { echo -e "\nPLS execute \". s.ssmuse.dot netcdff\"\n" ; false ; }
	s.locate --lib netcdff_s 1> /dev/null || { echo -e "\nPLS execute \". s.ssmuse.dot netcdff-4.4\"\n" ; false ; }

# RDIAG Diagnostic toolkit recipe

rdiag: initial_base
	echo "*** Making libdiag_sq98.a and libdiag_sq98_g.a ***" ;\
	cd $(DIAGNOSTIQUE)/src/lssub ; $(MAKE) VGDLIB=$(VGDLIB) ENTETE=$(ENTETE)
	echo "Making libprog_sq98.a" ;\
	cd $(DIAGNOSTIQUE)/src/lspgm ; $(MAKE)
	echo "*** Making executable r.diag ***" ;\
	cd $(DIAGNOSTIQUE)/src/lspgm ; $(MAKE) $(BASE_ARCH) OBJ="$(FIXES)" \
	NOPLOT=$(NOPLOT) GRAFLIB=$(GRAFLIB) DDFUN90=$(DDFUN90) VGDLIB=$(VGDLIB) \
	DIAG_VERSION=$(DIAG_VERSION) RMNLIB=$(RMNLIB) ENTETE=$(ENTETE)

# NetCDF to/from ( CCCma or CMC/RPN) file format converter recipe

cdf2conv: initial_base initial_cdf
	echo "*** Making libcdf2ccc.a ***" ;\
	cd $(DIAGNOSTIQUE)/src/cdf2ccc ; $(MAKE)
	echo "*** Making executable cdf2ccc ***" ;\
	cd $(DIAGNOSTIQUE)/src/cdf2ccc ;\
	$(MAKE) cdf2rpn CONV_VERSION=$(CONV_VERSION) \
	RMNLIB=$(RMNLIB) VGDLIB=$(VGDLIB) OBJ="$(FIXES)" \
	lNetCDF="$(lNetCDF)" UDUNITS="$(UDUNITS)" \
	DDFUN90=$(DDFUN90) ENTETE=$(ENTETE)
#	EXTRAS=$(EXTRAS)/NetcdfUdunits/$(EC_ARCH) DDFUN90=$(DDFUN90)

# Only generate the LSSUB, LSPM and CDF2CCC libraries

libs: initial_base initial_cdf
	echo "*** Making libdiag_sq98.a and libdiag_sq98_g.a ***" ;\
	cd $(DIAGNOSTIQUE)/src/lssub ; $(MAKE) VGDLIB=$(VGDLIB)
	echo "Making libprog_sq98.a" ;\
	cd $(DIAGNOSTIQUE)/src/lspgm ; $(MAKE) ;\
	echo "*** Making libcdf2ccc.a ***" ;\
	cd $(DIAGNOSTIQUE)/src/cdf2ccc ; $(MAKE)

# Online documentation (which was originaly found in $ARMNLIB/man/pdoc) recipe

document:
	cd $(DIAGNOSTIQUE)/man/pdoc ; $(MAKE) $@

web_document:
	cd $(DIAGNOSTIQUE)/man/pdoc ; $(MAKE) $@ \
	HOSTWEB=$(HOSTWEB) DIAGWEB=$(DIAGWEB)

# Clean

clean:
	cd $(DIAGNOSTIQUE)/src/lspgm   ; $(MAKE) $@ ;\
	cd $(DIAGNOSTIQUE)/src/lssub   ; $(MAKE) $@
