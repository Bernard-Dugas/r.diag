.SUFFIXES:

.SUFFIXES : .o .F90 .c .a

ifeq ($(BASE_ARCH),)
$(error FATAL: BASE_ARCH is not defined or empty, ABORTING)
endif

ifeq ($(BASE_ARCH),$(EC_ARCH))
	$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

SHELL = /bin/bash

ARUFLAG = U
AR      = s.ar
ARR     = $(AR) rv$(ARUFLAG)

LIB          = libdescrip.a

DIAGNOSTIQUE = $(CURDIR)/../../..
LIBDIR       = $(DIAGNOSTIQUE)/lib/$(EC_ARCH)
INCLUDEDIR   = $(DIAGNOSTIQUE)/include
MODULEDIR    = $(INCLUDEDIR)/$(EC_ARCH)

OBJECTS = $(EC_ARCH)/vgrid.o	$(EC_ARCH)/vgrid_utils.o	$(EC_ARCH)/vgrid_descriptors.o
MODFILES = $(EC_ARCH)/vgrid_utils.mod $(EC_ARCH)/vgrid_descriptors.mod

DEFINE = -DFICHSTD98

SUPP_OPT = -openmp

ifeq ($(USE_EFENCE),1)
   LEFENCE=-lefence
endif

include version.dot

%.o: %.F90
	s.f90 -c $(SUPP_OPT) -I$(EC_ARCH) -O0 -g $(MY_FFLAGS)  $(FFLAGS)$(LEFENCE) $<
	/bin/mv $*.o $*.mod $(EC_ARCH)

%.o: %.c
	s.cc -c $(SUPP_OPT) -I$(EC_ARCH) -O0 -g $(MY_CFLAGS) $(CFLAGS) $(LEFENCE) $<
	/bin/mv $*.o $(EC_ARCH)

default: updatelib

.PHONY:	updatelib	libdescrip	install	vgrid_version
.PHONY:	clean	distclean	clean_tilde	clean_all

dependencies.mk: *.c *.F90 *.h *.hc
	ls *.hc *.c *.F90 | ./s.dependencies.pl > dependencies.mk 

include dependencies.mk

updatelib:	libdescrip
	rsync -a $(EC_ARCH)/$(LIB) $(LIBDIR)
	rsync -a vgrid_version.h* $(EC_ARCH)/*.mod $(MODULEDIR)
	$(MAKE)	clean_all

libdescrip: install	vgrid_version $(OBJECTS) $(MODFILES)
	cd $(EC_ARCH) ; $(ARR) $(LIB) $(OBJECTS)

install:
	mkdir -p $(EC_ARCH)
	mkdir -p $(LIBDIR) 
	mkdir -p $(MODULEDIR)

vgrid.o:	vgrid.c \
	BODY_C_compute_heights_4001.hc  BODY_C_compute_pressure_2001.hc  armnlib.h  BODY_C_compute_pressure_5100.hc  \
	BODY_Cvgd_diag_withref_2ref.hc BODY_C_compute_pressure_5002_5003_5004_5005.hc  vgrid.h  \
	BODY_C_compute_heights_0001.hc  BODY_C_compute_pressure_1003_5001.hc  BODY_C_compute_pressure_1001_1002.hc   \
	BODY_C_compute_heights_21001.hc  

vgrid_descriptors.o:	vgrid_descriptors.F90 \
	vgrid_utils.o  BODY_F_dpidpis_withref.hf  BODY_F_dpidpis_withref_prof.hf  BODY_F_diag_withref.hf  \
	BODY_F_levels_withref_prof.hf   BODY_F_levels_withref.hf  BODY_F_diag_withref_prof.hf  vgrid_descriptors.hf  

vgrid_version:
	echo "character(len=128) :: vgrid_descriptors_version = \"VGRID $(BH_PULL_SOURCE_GIT_BRANCH) $(COMP_ARCH) $(ORDENV_PLAT) $(shell date)\"" > vgrid_version.hf
	echo                "char vgrid_descriptors_version[] = \"VGRID $(BH_PULL_SOURCE_GIT_BRANCH) $(COMP_ARCH) $(ORDENV_PLAT) $(shell date)\";" > vgrid_version.hc

vgrid_utils.o:	vgrid_utils.F90 vgrid_descriptors.hf  

clean:
	rm -f $(OBJECTS) $(MODFILES)

distclean: clean
	rm -rf $(EC_ARCH) ; rm -f vgrid_version.hc vgrid_version.hf dependencies.mk

clean_tilde:
	rm -f *~

clean_all: distclean  clean_tilde
