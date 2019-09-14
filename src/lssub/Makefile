.SUFFIXES :
.SUFFIXES : .cdk90 .ftn .ftn90 .txt .c .F90 .f .f90 .mod .o .a

# Description: Makefile principal de LSSUB

# Revision 3.71 2019/08/12 18:18  dugas
# - Utiliser libdescrip qui se trouve dans src/extras
#
# Revision 3.70 2019/05/22 16:16  dugas
# - Ajouter la variable ARUFLAG (=U, par defaut)
#
# Revision 3.69 2019/05/14 15:15  dugas
# - Ajouter le module averages_common utilise par le programme AVERAGES
#
# Revision 3.68 2019/02/19 13:13  dugas
# - Enlever toutes re-definitions de MAKE qui pourraient inclure
#   une revision de la valeur de la variable DIAGNOSTIQUE. Cette
#   variable est deja correctement definie localement
#
# Revision 3.67 2019/02/12 14:14  dugas
# - Remplacer les instances de MAKE par des MYMAKE
#
# revision 3.66 2017/11/02 16:16  dugas
# - Modifier les recettes associees aux modules
# - Ajouter un segment "make module" dans la cible "initial"
# - Modifier/Completer les dependances aux modules de r.diag
#
# Revision 3.65  2016/10/27 13:24  dugas
# - Le macro ENTETE passe de 28 a 32
# - Modifier la recette .cdk.o ("rm" plutot que "mv" le .f90)
# - Ajouter le traitement de diag_convip_plus.ftn90
# - Corriger la syntaxe de la liste de dependance de convsub.o et la midifier
# - Modifier les dependances de coord.o, diag_convert_ip123.o, diag_toc.o, get.o et util.o
# - Modifier le second make de la cible principale (on definit le macro OPTIM90)

ifeq ($(BASE_ARCH),)
$(error FATAL: BASE_ARCH is not defined or empty, ABORTING)
endif

ifeq "$(BASE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

SHELL = /bin/bash

ARUFLAG = U
AR      = s.ar
ARR     = $(AR) rv$(ARUFLAG)

CPPT    = r.gppf -chop_bang -gpp -F
ZIP     = $(shell which zip)

LSSDIR  = $(CURDIR)
DIAGNOSTIQUE = $(LSSDIR)/../..
MYMAKE  = $(MAKE) DIAGNOSTIQUE=$(DIAGNOSTIQUE) ARUFLAG=$(ARUFLAG)

CFLAGS  =
FFLAGS  =

HOSTWEB = 
DIAGWEB =

ENTETE  = 32
NP      = 1
OPTIL   = 2
STD     = 98

DIAGS   = $(DIAGNOSTIQUE)
LIBDIR  = $(DIAGS)/lib/$(EC_ARCH)
MANDIR  = $(DIAGS)/man/pdoc
INCLUDE = $(DIAGS)/include
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

.cdk90.o:
	$(FCOMP) $<
	/bin/rm -f $*.f90

#.cdk90.o:
#	$(FCOMP) $<
#	/bin/mv -f $*.f90 $(SOURCES)

.c.o:
	$(CC) $<

.cdk90.mod:
	$(FTNC) $<
	$(FCOMP) $*.f90
	rsync -a $@ $(MODDIR)
	/bin/rm -f $@ $*.f90 $*.o

.f90.mod:
	$(FCOMP) $<
	rsync -a $@ $(MODDIR)
	/bin/rm -f $@ $*.o

.ftn.f:
	$(FTNC) $<

.cdk90.f90:
	$(FTNC) $<

.ftn90.f90:
	$(FTNC) $<

.ftn.txt:
	$(FTNT) $< | sed s/'^[*]'/' '/ > $(MANDIR)/$@ ;\
	if ( test "$(HOSTWEB):$(DIAGWEB)" != ":" ) then \
	rsync -av $(MANDIR)/$@ $(HOSTWEB):$(DIAGWEB) ; fi

.f.a:
	$(FCOMP) $<
	if [ -f $@ ]; then /bin/chmod 644 $@ ; fi
	$(ARR) $@ $*.o
	/bin/cp -p $*.f $(SOURCES)
	/bin/rm -f $*.o

.F90.o:
	s.f90 -c $<

.f90.a:
	$(FCOMP) $<
	if [ -f $@ ]; then /bin/chmod 644 $@ ; fi
	$(ARR) $@ $*.o
	/bin/cp -p $*.f90 $(SOURCES)
	/bin/rm -f $*.o

.o.a:
	$(ARR) $@ $*.o

%.o: %.mod

# Listes des fichiers (sources, documentation, objets)

OBJECTS= \
	alp.o		cofg.o		convsub.o 	coord.o		ddl.o \
	diag_convert_ip123.o	diag_convip_plus.o	diag_toc.o \
	dimgt.o	divers.o 	eapl.o		eof.o		fastaf.o \
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
	gaussg.F90	util2.F90

f90DECKS= \
	diag_convert_ip123.f90	diag_convip_plus.f90 \
	diag_toc.f90	filter.f90	averages_common.f90 \
	stats_signatures.f90

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
	filter.ftn90	diag_convip_plus.ftn90

CDK90DECKS = \
	diag_convert_ip123.cdk90	diag_toc.cdk90 \
	stats_signatures.cdk90		averages_common.cdk90

MODDECKS= \
	$(MODDIR)/diag_convert_ip123.mod 	$(MODDIR)/diag_toc.mod \
	$(MODDIR)/stats_signatures.mod		$(MODDIR)/averages_common.mod

# Dependances

alp.f	: alp.ftn \
	$(INCLUDE)/calpi.cdk
convsub.f : convsub.ftn \
	$(INCLUDE)/fstparm.cdk $(INCLUDE)/highbuf.cdk \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk \
	$(MODDIR)/diag_convert_ip123.mod \
	$(MODDIR)/diag_toc.mod
coord.f   : coord.ftn \
	$(INCLUDE)/highbuf.cdk
fastaf2.f : fastaf2.ftn \
	$(INCLUDE)/diagmacros.cdk
get.f     : get.ftn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk \
	$(MODDIR)/stats_signatures.mod
inctdia.f : inctdia.ftn \
	$(INCLUDE)/consdia.cdk
label.f   : label.ftn \
	$(INCLUDE)/machtype.cdk \
	$(MODDIR)/diag_convert_ip123.mod \
	$(MODDIR)/stats_signatures.mod
lgtst.f	: lgtxt.ftn \
	$(INCLUDE)/calpi.cdk
lowio92.f : lowio92.ftn \
	$(INCLUDE)/machtype.cdk
pacc92.f  : pacc92.ftn \
	$(INCLUDE)/machtype.cdk
paccrn.f  : paccrn.ftn \
	$(INCLUDE)/machtype.cdk
put.f     : put.ftn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/highbuf.cdk \
	$(MODDIR)/diag_convert_ip123.mod \
	$(MODDIR)/stats_signatures.mod
rec.f     : rec.ftn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk
util.f    : util.ftn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk

averages_common.f90	: averages_common.cdk90 \
	$(MODDIR)/stats_signatures.mod
diag_convip_plus.f90 : diag_convip_plus.ftn90 \
	$(MODDIR)/diag_convert_ip123.mod

averages_common.mod	: averages_common.cdk90 \
	$(MODDIR)/stats_signatures.mod

$(LIBS)(averages_common.o) : averages_common.cdk90 \
	$(MODDIR)/stats_signatures.mod
$(LIBS)(convsub.o)  : convsub.ftn \
	$(INCLUDE)/fstparm.cdk	$(INCLUDE)/highbuf.cdk \
	$(INCLUDE)/machtype.cdk	$(INCLUDE)/zlabel.cdk \
	$(MODDIR)/diag_convert_ip123.mod \
	$(MODDIR)/diag_toc.mod
$(LIBS)(diag_convip_plus.o) : diag_convip_plus.ftn90 \
	$(MODDIR)/diag_convert_ip123.mod
$(LIBS)(get.o)      : get.ftn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk \
	$(MODDIR)/diag_convert_ip123.mod \
	$(MODDIR)/stats_signatures.mod
$(LIBS)(label.o)    : label.ftn \
	$(MODDIR)/diag_convert_ip123.mod \
	$(MODDIR)/stats_signatures.mod
$(LIBS)(put.o)      : put.ftn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk \
	$(MODDIR)/diag_convert_ip123.mod \
	$(MODDIR)/stats_signatures.mod
$(LIBS)(util.o)     : util.ftn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk

# Liste des cibles

initial:
# Setup des repertoires et liens dependants de $(EC_ARCH)
	/bin/mkdir -p $(SOURCES) $(LIBDIR) $(MANDIR) $(MODDIR)
#	s.locate --lib $(VGDLIB) 1> /dev/null || { echo -e \n PLS execute \". s.ssmuse.dot vgriddesc\" \n ; false ; }
	if [[ ! -f $(LIBDIR)/libdescrip.a || -z "$(VGDLIB)" ]]; then \
		cd $(DIAGNOSTIQUE)/src/extras/descrip ; $(MAKE) ARUFLAG=$(ARUFLAG) ; fi
	if [[ `/bin/ls -L $(SOURCES)/Makefile ; echo $?` != 0 ]]; then \
	/bin/ln -sf $(DIAGS)/src/lssub/Makefile_sources$(STD) $(SOURCES)/Makefile ; fi
	if [[ `/bin/ls -L $(MODDIR)/Makefile ; echo $?` != 0 ]]; then \
	/bin/ln -sf $(INCLUDE)/Makefile_mods $(MODDIR)/Makefile ; fi

# Voir "Managing project with make" de O'Reilly & Associates, 1989,
# pages 51-53, pour les details de ce qui se passe dans ce qui suit
$(LIBS):	initial	$(FTNDECKS)	$(FTN90DECKS) \
	$(F90DECKS)	$(CDECKS)	$(CDK90DECKS) \
	$(MODDECKS) \
	$(LIBS)(alp.o)		\
	$(LIBS)(averages_common.o)		\
	$(LIBS)(cofg.o)		\
	$(LIBS)(convsub.o)	\
	$(LIBS)(coord.o)		\
	$(LIBS)(ddl.o)		\
	$(LIBS)(diag_convert_ip123.o)		\
	$(LIBS)(diag_convip_plus.o)		\
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
	$(LIBS)(stats_signatures.o)	\
	$(LIBS)(temperton.o) 	\
	$(LIBS)(trans1d.o) 	\
	$(LIBS)(util.o)		\
	$(LIBS)(util2.o)		\
	$(LIBS)(vfft.o)

clean:
# Faire un petit menage. On enleve tous les .o
	/bin/rm -f *.o *.mod

$(ARCH): initial
# Make en trois etapes: 1) Le premier $(MAKE) mets a jour les modules.
# 2) Les deux suivant $(MAKE) mettent a jour les sources dans $(SOURCES),
# tandis que 3) les deux derniers generent les librairies (ver. debug et
# et ver. ordinaire, respectivement)
	cd $(MODDIR)  ; $(MAKE) ;\
	cd $(LSSDIR)  ; $(MAKE) $(LIBS) ; $(MAKE) clean ;\
	cd $(SOURCES) ; $(MYMAKE) ; $(MYMAKE) DIAGLIB=diag_sq$(STD) OPTIM='-O 2' OPTIM90='-O 2'
