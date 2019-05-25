
ifeq ($(BASE_ARCH),)
$(error FATAL: BASE_ARCH is not defined or empty, ABORTING)
endif

ifeq "$(BASE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

DIAGNOSTIQUE = $(CURDIR)/../../..
LIBDIR       = $(DIAGNOSTIQUE)/lib/$(EC_ARCH)
INCLUDEDIR   = $(DIAGNOSTIQUE)/include/$(EC_ARCH)
LIB          = libudunits2f_s.a

default: lib

.PHONY: lib

lib: install f_udunits_2.f90 f_udunits_2.inc
	s.f90 -c -O f_udunits_2.f90
	s.ar rcv $(LIBDIR)/$(LIB) f_udunits_2.o
	/bin/mv *.mod $(INCLUDEDIR)
	/bin/rm f_udunits_2.o

install:
	mkdir -p $(LIBDIR)
	mkdir -p $(INCLUDEDIR)
