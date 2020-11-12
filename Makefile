
# RDIAG Toolkit main Makefile

ifeq ($(BASE_ARCH),)
$(error FATAL: BASE_ARCH is not defined or empty, ABORTING)
endif

ifeq "$(BASE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

SHELL   = /bin/bash

# Location of the development versions of the src, bin/$(BASE_ARCH),
# lib/$(EC_ARCH) and man/pdoc directory trees

DIAGNOSTIQUE = $(CURDIR)

# Destination that the current working binaries and libraries
# will be exported to, as well as where to find the EXTRAS

DESTINATION = $(DIAGNOSTIQUE)/..

# NetCDF v_3.6 and Udunits v_1.2 libraries used to be found 
# in the EXTRAS directory. The package now can either 1) use
# the netcdff SSM package to locate all the necessary shared
# libraries including MFV's udunits2f FORTRAN wrapper as well
# as udunits2 itself or 2) use the static versions of these
# same libraries as provided by the netcdff-4.4 SSM package.
# As of May 2017, the second option is always used. And as
# of May 2019, the udunits2f package is included in these
# extras.
#EXTRAS  = $(DESTINATION)/extras
#EXTLIB  = $(EXTRAS)/NetcdfUdunits/$(EC_ARCH)/lib

# Directories used/created by this Makefile

BINDIR  = $(DIAGNOSTIQUE)/bin/$(BASE_ARCH)
INCLUDE = $(DIAGNOSTIQUE)/include
LIBDIR  = $(DIAGNOSTIQUE)/lib/$(EC_ARCH)
MANDIR  = $(DIAGNOSTIQUE)/man/pdoc
MODDIR  = $(INCLUDE)/$(EC_ARCH)
SUBDIR  = $(DIAGNOSTIQUE)/src/lssub/sources$(STD)/$(EC_ARCH)

# Doit-on ajouter le parametre U a la commande s.ar ? Par defaut, oui !

ARUFLAG = U

# Include (very) old NCAR graphics (default=no) ?

NOPLOT  = -DNOPLOT
GRAFLIB =

# Binaires pre-compiles

#FIXES   = $(PWD)/lib/$(EC_ARCH)/stubs.o

# WEB Host Server for the documentation

HOSTWEB = pascal
DIAGWEB = public_html

# RMN and Vgrid_Descriptor library names

RMNLIB  = rmn_016.3.1
VGDLIB  = descrip

# LAPCK and BLAS libraries

BLAS    = blas
LAPACK  = lapack

# DDFUN90 (Substitute Quad-Precision) library name

DDFUN90  = ddfun90

ifeq ($(SHARED_NETCDF),)
# Static (_s) load via symlinks in the netcdff-4.4 SSM package
NLocate  = s.locate
lNetCDF  = netcdff_s netcdf_s hdf5_hl_s hdf5_s dl sz_s z
UDUNITS  = udunits2f_s udunits2_s expat
else
# Dynamic (shared-object) load via the system's netcdff package
NLocate  = true
lNetCDF  =
UDUNITS  = udunits2f_s udunits2 expat
endif

DIAG_VERSION = 6.4.5
CONV_VERSION = 2.3.5

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
	if [[ `/bin/ls -L $(MODDIR)/Makefile ; echo $?` != 0 ]]; then \
		/bin/ln -sf $(INCLUDE)/Makefile_mods $(MODDIR)/Makefile ; fi
#	s.locate --lib $(VGDLIB) 1> /dev/null || { echo -e \n PLS execute \". s.ssmuse.dot vgriddesc\" \n ; false ; }
	if [[ ! -f $(LIBDIR)/libdescrip.a || -z "$(VGDLIB)" ]]; then \
		cd $(DIAGNOSTIQUE)/src/extras/descrip ; $(MAKE) ARUFLAG=$(ARUFLAG) ; fi
	if [[ ! -f $(LIBDIR)/libddfun90.a || -z "$(DDFUN90)" ]]; then \
		cd $(DIAGNOSTIQUE)/src/extras/ddfun90 ; $(MAKE) RMNLIB=$(RMNLIB) ; fi
	if [[ ! -f $(LIBDIR)/program_version.o ]]; then cd $(LIBDIR) ;	s.f77 -g -c ../../program_version.f ; fi
	if [[ ! -x $(BINDIR)/r.echo ]]; then cd $(DIAGNOSTIQUE)/src/extras/tools ; $(MAKE) ; fi
	if [[ `/usr/bin/diff $(DIAGNOSTIQUE)/bin/r.diag_commands \
		$(BINDIR)/r.diag_commands 1> /dev/null 2>&1 ; echo $?` != 0 ]]; then \
		rsync -a $(DIAGNOSTIQUE)/bin/r.diag_commands $(BINDIR) ; fi

initial_cdf:
	$(NLocate) --lib netcdff_s 1> /dev/null || { echo -e "\nPLS execute \". s.ssmuse.dot netcdff-4.4\"\n" ; false ; }
	if [[ ! -f $(LIBDIR)/libudunits2f_s.a ]]; then cd $(DIAGNOSTIQUE)/src/extras/udunits-f-2.0 ; $(MAKE) ; fi

# RDIAG Diagnostic toolkit recipe

rdiag: initial_base
	echo "*** Making the RDIAG modules ***" ; cd $(MODDIR) ; $(MAKE)
	echo "*** Making libdiag_sq98.a and libdiag_sq98_g.a ***" ;\
	cd $(DIAGNOSTIQUE)/src/lssub ; $(MAKE) VGDLIB=$(VGDLIB) ENTETE=$(ENTETE) ARUFLAG=$(ARUFLAG)
	echo "Making libprog_sq98.a" ;\
	cd $(DIAGNOSTIQUE)/src/lspgm ; $(MAKE) ARUFLAG=$(ARUFLAG)
	echo "*** Making executable r.diag ***" ;\
	cd $(DIAGNOSTIQUE)/src/lspgm ; $(MAKE) $(BASE_ARCH) OBJ="$(FIXES)" \
	NOPLOT=$(NOPLOT) GRAFLIB=$(GRAFLIB) DDFUN90=$(DDFUN90) VGDLIB=$(VGDLIB) \
	DIAG_VERSION=$(DIAG_VERSION) RMNLIB=$(RMNLIB) ENTETE=$(ENTETE) \
	BLAS=$(BLAS) LAPACK=$(LAPACK)

# NetCDF to/from ( CCCma or CMC/RPN) file format converter recipe

cdf2conv: initial_base initial_cdf
	echo "*** Making libcdf2ccc.a ***" ;\
	cd $(DIAGNOSTIQUE)/src/cdf2ccc ; $(MAKE) ARUFLAG=$(ARUFLAG) NLocate=$(NLocate)
	echo "*** Making executable cdf2ccc ***" ;\
	cd $(DIAGNOSTIQUE)/src/cdf2ccc ;\
	$(MAKE) cdf2rpn CONV_VERSION=$(CONV_VERSION) \
	RMNLIB=$(RMNLIB) VGDLIB=$(VGDLIB) OBJ="$(FIXES)" \
	lNetCDF="$(lNetCDF)" UDUNITS="$(UDUNITS)" \
	DDFUN90=$(DDFUN90) ENTETE=$(ENTETE) NLocate=$(NLocate)

# Only generate the LSSUB, LSPM and CDF2CCC libraries

libs: initial_base initial_cdf
	echo "*** Making libdiag_sq98.a and libdiag_sq98_g.a ***" ;\
	cd $(DIAGNOSTIQUE)/src/lssub ; $(MAKE) VGDLIB=$(VGDLIB) ENTETE=$(ENTETE) ARUFLAG=$(ARUFLAG)
	echo "Making libprog_sq98.a" ;\
	cd $(DIAGNOSTIQUE)/src/lspgm ; $(MAKE) ARUFLAG=$(ARUFLAG) ;\
	echo "*** Making libcdf2ccc.a ***" ;\
	cd $(DIAGNOSTIQUE)/src/cdf2ccc ; $(MAKE) ARUFLAG=$(ARUFLAG)

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

veryclean: clean
	/bin/rm -rf $(SUBDIR) $(MODDIR) $(LIBDIR) $(BINDIR)
