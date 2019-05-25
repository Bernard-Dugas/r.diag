
ifeq ($(BASE_ARCH),)
$(error FATAL: BASE_ARCH is not defined or empty, ABORTING)
endif

ifeq "$(BASE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

DIAGNOSTIQUE = $(CURDIR)/../../..

F_SRC = ddfun90.f90  ddmod90.f90

F_OBJ = ddfun90.o  ddmod90.o

C_SRC = x86.c

C_OBJ = x86.o

RMNLIB = rmn_015

LIB = libddfun90.a

INCDIR = $(DIAGNOSTIQUE)/include/$(EC_ARCH)
LIBDIR = $(DIAGNOSTIQUE)/lib/$(EC_ARCH)

default: lib

.PHONY: lib

lib:	$(LIBDIR)/$(LIB)

$(LIBDIR)/$(LIB):
	mkdir -p $(INCDIR) $(LIBDIR)
	rm -f *.o *.mod
	s.f90 -c -O $(F_SRC)
	s.cc -c $(C_SRC)
	ar rcv $@ *.o
	rm -f *.o
	mv *.mod $(INCDIR)

test:	$(LIBDIR)/$(LIB)
	rm -f a.out
	s.f90 -O tquadtsdd.f90 -I$(INCDIR) -L$(LIBDIR) -lddfun90 -l$(RMNLIB)
	./a.out
	rm a.out

clean:
	rm -f *.o *.mod a.out
	rm -rf $(LIBDIR)/* $(INCDIR)/*
