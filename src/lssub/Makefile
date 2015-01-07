.SUFFIXES :
.SUFFIXES : .cdk90 .ftn .ftn90 .txt .c .F90 .f .f90 .mod .o .a

# LSSUB Makefile

ifeq "$(BASE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

SHELL = /bin/bash

AR      = s.ar
CPPT    = r.gppf -chop_bang -gpp -F
ZIP     = `which zip`

DIAGNOSTIQUE = $(CURDIR)/../..
MAKE    = make DIAGNOSTIQUE=$(DIAGNOSTIQUE)

CFLAGS  =
FFLAGS  = 

ENTETE  = 28
NP      = 1
OPTIL   = 2
STD     = 98

DIAGS   = $(DIAGNOSTIQUE)
INCLUDE = $(DIAGS)/include
LIBDIR  = $(DIAGS)/lib/$(EC_ARCH)
MANDIR  = $(DIAGS)/man/pdoc
MODDIR  = $(INCLUDE)/$(EC_ARCH)
SOURCES = $(DIAGS)/src/lssub/sources$(STD)/$(EC_ARCH)

DIAGLIB = diag_sq$(STD)
LIBS    = $(LIBDIR)/lib$(DIAGLIB).a
VGDLIB  = descrip

TXT     = -DNO_SOURCE -DAUTO_DOC
DEFINE  = -Dtaille_entete=$(ENTETE) -Dnombre_de_taches=$(NP)

COMPILE = s.compile
FCOMP   = $(COMPILE) -includes "$(INCLUDE) $(MODDIR)" -defines "=$(DEFINE)" -O $(OPTIL) -optf "$(FFLAGS)" -src
CC      = $(COMPILE) -includes "$(INCLUDE)"           -defines "=$(DEFINE)" -O $(OPTIL) -optc "$(CFLAGS)" -src
FTNC    = $(COMPILE) -includes "$(INCLUDE)"           -defines "=$(DEFINE)" -P -src
FTNT    = $(CPPT) $(TXT)

ARCH    = $(BASE_ARCH)

default: $(ARCH)

.PRECIOUS: $(LIBS)

.DEFAULT:
	cvs checkout $@

.cdk90.a:
	/bin/rm -f $*.f90
	$(FTNC) $<
	/bin/mkdir -p $(SOURCES)
	/bin/mv -f $*.f90 $(SOURCES)

.F90.a:
	/bin/mkdir -p $(SOURCES)
	rsync -a  $< $(SOURCES)

.ftn.a:
	/bin/rm -f $*.f
	$(FTNC) $<
	/bin/mv -f $*.f $(SOURCES)

.ftn90.a:
	/bin/rm -f $*.f90
	$(FTNC) $<
	/bin/mkdir -p $(SOURCES)
	/bin/mv -f $*.f90 $(SOURCES)

.c.a:
	/bin/cp -p $< $(SOURCES)

.ftn.o:
	$(FCOMP) $<
	/bin/rm -f $*.f

.ftn90.o:
	$(FCOMP) $<
	/bin/rm -f $*.f90

.c.o:
	$(CC) $<

.cdk90.mod:
	$(FCOMP) $<
	rsync -a $@ $(MODDIR)

.f90.mod:
	$(FCOMP) $<
	rsync -a $@ $(MODDIR)

.cdk90.o:
	$(FCOMP) $<
	/bin/mv -f $*.f90 $(SOURCES)

.ftn.f:
	$(FTNC) $<

.cdk90.f90:
	$(FTNC) $<

.ftn90.f90:
	$(FTNC) $<

.ftn.txt:
	/bin/rm -f $(MANDIR)/$*.txt $*.txt
	$(FTNT) $< | sed s/'^[*]'/' '/ > $(MANDIR)/$*.txt
	rsync -av $(MANDIR)/$*.txt .
	if ( test "$(HOSTWEB):$(DIAGWEB)" != ":" ) then \
	rsync -av $(MANDIR)/$*.txt $(HOSTWEB):$(DIAGWEB) ; fi

.f.a:
	$(FCOMP) $<
	if [ -f $@ ]; then /bin/chmod 644 $@ ; fi
	$(AR) rv $@ $*.o
	/bin/cp -p $*.f $(SOURCES)
	/bin/rm -f $*.o

.F90.o:
	s.f90 -c $<

.f90.a:
	$(FCOMP) $<
	if [ -f $@ ]; then /bin/chmod 644 $@ ; fi
	$(AR) rv $@ $*.o
	/bin/cp -p $*.f90 $(SOURCES)
	/bin/rm -f $*.o

.o.a:
	$(AR) rv $@ $*.o

%.o: %.mod

# Listes des fichiers (sources, documentation, objets)

OBJECTS= \
	alp.o		cofg.o		convsub.o 	coord.o		ddl.o \
	dimgt.o	diag_toc.o	 	divers.o 	eapl.o		eof.o		fastaf.o \
	fastaf2.o 	ffgfw2.o 	fftcray.o 	fftini.o		ffwfg2.o \
	filev.o	 	filter.o 	four2.o		gaussg.o		gcof.o \
	gemdiag.o 	get.o		ggd.o		gobits.o 	inctdia.o \
	interp.o		ism.o		label.o		lgtst.o		lire_arg.o \
	lowio.o	 	lowio92.o	lssub.o		lwbw.o		mem.o \
	mrcdiag.o	newtim.o  	paccrn.o 	pacc92.o		pael.o \
	pfa.o		precon.o	prob.o	prob.o		put.o		rec.o \
	sfft.o		spectra.o 	temperton.o	trans1d.o	util.o \
	util2.o	vfft.o

CDECKS= \
	divers.c

HTMLDECKS = Error_messages.html

FDECKS= \
	alp.f		cofg.f		convsub.f 	coord.f		ddl.f \
	dimgt.f		eapl.f		eof.f		fastaf.f \
	fastaf2.f 	ffgfw2.f 	fftcray.f 	fftini.f		ffwfg2.f \
	filev.f	 	four2.f		gcof.f \
	gemdiag.f 	get.f		ggd.f		gobits.f		inctdia.f \
	interp.f		ism.f		label.f		lgtst.f		lire_arg.f \
	lowio.f		lowio92.f	lssub.f		lwbw.f		mem.f \
	mrcdiag.f	newtim.f  	paccrn.f 	pacc92.f		pael.f \
	pfa.f		precon.f	prob.f 	prob.o		put.f		rec.f \
	sfft.f		spectra.f 	temperton.f	trans1d.f	util.f \
	vfft.f

F90DECKS= \
	util2.F90

f90DECKS= \
	diag_convert_ip123.f90 	diag_toc.f90	filter.f90	gaussg.f90

FTNDECKS= \
	alp.ftn		cofg.ftn		convsub.ftn 	coord.ftn 	ddl.ftn \
	dimgt.ftn 			eapl.ftn		eof.ftn		fastaf.ftn \
	fastaf2.ftn 	ffgfw2.ftn	fftcray.ftn 	fftini.ftn	ffwfg2.ftn \
	filev.ftn 	four2.ftn	gcof.ftn \
	gemdiag.ftn	get.ftn		ggd.ftn		gobits.ftn	inctdia.ftn \
	interp.ftn	ism.ftn		label.ftn 	lgtst.ftn 	lire_arg.ftn \
	lowio.ftn	lowio92.ftn	lssub.ftn 	lwbw.ftn 	mem.ftn \
	mrcdiag.ftn	newtim.ftn  	paccrn.ftn 	pacc92.ftn	pael.ftn \
	pfa.ftn		precon.ftn	prob.ftn 	put.ftn  	rec.ftn \
	sfft.ftn		spectra.ftn 	temperton.ftn	trans1d.ftn	util.ftn \
	vfft.ftn

FTN90DECKS = \
	filter.ftn90	gaussg.ftn90

CDK90DECKS = \
	diag_convert_ip123.cdk90	diag_toc.cdk90

MODDECKS= \
	diag_convert_ip123.mod 	diag_toc.mod

# Dependances

alp.f	: alp.ftn \
	$(INCLUDE)/calpi.cdk
convsub.f : convsub.ftn \
	$(INCLUDE)/fstparm.cdk $(INCLUDE)/highbuf.cdk \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk
coord.f   : coord.ftn \
	$(INCLUDE)/highbuf.cdk
fastaf2.f : fastaf2.ftn \
	$(INCLUDE)/diagmacros.cdk
get.f     : get.ftn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk
inctdia.f : inctdia.ftn \
	$(INCLUDE)/consdia.cdk
label.f   : label.ftn \
	$(INCLUDE)/machtype.cdk
lgtst.f	: lgtxt.ftn \
	$(INCLUDE)/calpi.cdk
lowio92.f : lowio92.ftn \
	$(INCLUDE)/machtype.cdk
pacc92.f  : pacc92.ftn \
	$(INCLUDE)/machtype.cdk
paccrn.f  : paccrn.ftn \
	$(INCLUDE)/machtype.cdk
put.f     : put.ftn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/highbuf.cdk
rec.f     : rec.ftn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk
util.f    : util.ftn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk

convsub.o  : convsub.ftn diag_toc.mod \
	$(INCLUDE)/fstparm.cdk	$(INCLUDE)/highbuf.cdk \
	$(INCLUDE)/machtype.cdk	$(INCLUDE)/zlabel.cdk
coord.o : coord.ftn	diag_convert_ip123.mod
diag_convert_ip123.o : diag_convert_ip123.cdk90
diag_toc.o : diag_toc.cdk90	vgrid_descriptors.mod
get.o      : get.ftn diag_toc.mod \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk
util.o     : util.ftn  diag_toc.mod \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk

# Liste des cibles

initial:
# Setup des repertoires et liens dependants de $(EC_ARCH)
	/bin/mkdir -p $(SOURCES) $(LIBDIR) $(MANDIR) $(MODDIR)
	s.locate --lib $(VGDLIB) 1> /dev/null || { echo -e \n PLS execute \". s.ssmuse.dot vgriddesc\" \n ; false ; }
	if [[ `/bin/ls -L $(SOURCES)/Makefile ; echo $?` != 0 ]]; then \
	/bin/ln -sf $(DIAGS)/src/lssub/Makefile_sources$(STD) $(SOURCES)/Makefile ; fi

document: lssub.ftn
# Creation de la documentation pour lssub
	if ( test "$(HOSTWEB):$(DIAGWEB)" != ":" ) then \
	rsync -av $(HTMLDECKS) $(HOSTWEB):$(DIAGWEB) ; fi
	$(MAKE) lssub.txt

# Voir "Managing project with make" de O'Reilly & Associates, 1989,
# pages 51-53, pour les details de ce qui se passe dans ce qui suit
$(LIBS):	initial	$(FTNDECKS)	$(FTN90DECKS) \
	$(F90DECKS)	$(CDECKS)	$(CDK90DECKS) \
	$(MODDECKS) \
	$(LIBS)(alp.o)		\
	$(LIBS)(cofg.o)		\
	$(LIBS)(convsub.o)	\
	$(LIBS)(coord.o)		\
	$(LIBS)(ddl.o)		\
	$(LIBS)(diag_convert_ip123.o)		\
	$(LIBS)(diag_toc.o)		\
	$(LIBS)(dimgt.o)		\
	$(LIBS)(divers.o)		\
	$(LIBS)(eapl.o)		\
	$(LIBS)(eof.o)		\
	$(LIBS)(fastaf.o)		\
	$(LIBS)(fastaf2.o)	\
	$(LIBS)(ffgfw2.o)		\
	$(LIBS)(fftcray.o)	\
	$(LIBS)(fftini.o)		\
	$(LIBS)(ffwfg2.o)		\
	$(LIBS)(filev.o)		\
	$(LIBS)(filter.o)		\
	$(LIBS)(four2.o)		\
	$(LIBS)(gaussg.o)		\
	$(LIBS)(gcof.o)		\
	$(LIBS)(gemdiag.o)	\
	$(LIBS)(get.o)		\
	$(LIBS)(ggd.o)		\
	$(LIBS)(gobits.o)		\
	$(LIBS)(inctdia.o)	\
	$(LIBS)(interp.o)		\
	$(LIBS)(ism.o)		\
	$(LIBS)(label.o)		\
	$(LIBS)(lgtst.o)		\
	$(LIBS)(lire_arg.o)	\
	$(LIBS)(lowio.o)		\
	$(LIBS)(lowio92.o)	\
	$(LIBS)(lssub.o)		\
	$(LIBS)(lwbw.o)		\
	$(LIBS)(mem.o)		\
	$(LIBS)(mrcdiag.o)	\
	$(LIBS)(newtim.o)		\
	$(LIBS)(paccrn.o)		\
	$(LIBS)(pacc92.o)		\
	$(LIBS)(pael.o)		\
	$(LIBS)(pfa.o)		\
	$(LIBS)(precon.o)		\
	$(LIBS)(prob.o)		\
	$(LIBS)(put.o)		\
	$(LIBS)(rec.o) 		\
	$(LIBS)(sfft.o)		\
	$(LIBS)(spectra.o)	\
	$(LIBS)(temperton.o) 	\
	$(LIBS)(trans1d.o) 	\
	$(LIBS)(util.o)		\
	$(LIBS)(util2.o)		\
	$(LIBS)(vfft.o)

clean:
# Faire un petit menage. On enleve tous les .o
	/bin/rm -f *.o *.mod

$(ARCH): initial
# Make en deux etapes: Les deux premiers $(MAKE) mettent a jour les sources
# dans $(SOURCES) tandis que les deux derniers generent les librairies 
# (debug et ordinaire, respectivement)
	$(MAKE) $(LIBS) ;	$(MAKE) clean ; cd $(SOURCES) ; \
	$(MAKE) ; $(MAKE) DIAGLIB=diag_sq$(STD) OPTIM='-O 2'
