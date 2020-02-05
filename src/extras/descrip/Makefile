
ifeq ($(BASE_ARCH),)
$(error FATAL: BASE_ARCH is not defined or empty, ABORTING)
endif

ifeq ($(BASE_ARCH),$(EC_ARCH))
	$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

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

ifeq ($(USE_EFENCE),1)
   LEFENCE=-lefence
endif

include version.dot

%.o: %.F90
	s.f90 -c -openmp -I$(EC_ARCH) -O0 -g $(MY_FFLAGS)  $(FFLAGS)$(LEFENCE) $<
	/bin/mv $*.o $*.mod $(EC_ARCH)

%.o: %.c
	s.cc -c -openmp -I$(EC_ARCH) -O0 -g $(MY_CFLAGS) $(CFLAGS) $(LEFENCE) $<
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

vgrid_version:
	echo "character(len=128) :: vgrid_descriptors_version = \"VGRID $(BH_PULL_SOURCE_GIT_BRANCH) $(COMP_ARCH) $(ORDENV_PLAT) $(shell date)\"" > vgrid_version.hf
	echo                "char vgrid_descriptors_version[] = \"VGRID $(BH_PULL_SOURCE_GIT_BRANCH) $(COMP_ARCH) $(ORDENV_PLAT) $(shell date)\";" > vgrid_version.hc

clean:
	rm -f $(OBJECTS) $(MODFILES)

distclean: clean
	rm -rf $(EC_ARCH) ; rm -f vgrid_version.hc vgrid_version.hf dependencies.mk

clean_tilde:
	rm -f *~

clean_all: distclean  clean_tilde
