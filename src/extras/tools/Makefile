
ifeq ($(BASE_ARCH),)
$(error FATAL: BASE_ARCH is not defined or empty, ABORTING)
endif

DIAGNOSTIQUE = $(CURDIR)/../../..
BINDIR       = $(DIAGNOSTIQUE)/bin/$(BASE_ARCH)
CC           = s.cc

default: r.echo

r.echo:	r.echo.c
	$(CC) -o $@ $@.c
	/bin/mv $@ $(DIAGNOSTIQUE)/bin/$(BASE_ARCH)
	/bin/rm -f $@.o

clean:
	/bin/rm -f *.o 
