.SUFFIXES :
.SUFFIXES : .txt .ftn .f .h .o .a

# CDF2CCC Makefile

ifeq "$(BASE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

SHELL    = /bin/bash

AR       = s.ar
CPPT     = r.gppf -chop_bang -gpp -F

FFLAGS   =  
CFLAGS   =

FCOMPF   =  -debug -O 0
COMPF    = 

DIAGNOSTIQUE  = $(CURDIR)/../..
IGNORE_ERRORS = set -e
MAKE     = $(IGNORE_ERRORS) ; make DIAGNOSTIQUE=$(DIAGNOSTIQUE)

EXTRAS   = $(DIAGNOSTIQUE)/extras/NetcdfUdunits/$(EC_ARCH)

DIAGS    = $(DIAGNOSTIQUE)
BINDIR   = $(DIAGS)/bin/$(BASE_ARCH)
INCLUDES = $(DIAGS)/include
LIBDIR   = $(DIAGS)/lib/$(EC_ARCH)
MANDIR   = $(DIAGS)/man/pdoc
MODDIR   = $(INCLUDES) $(INCLUDES)/$(EC_ARCH) $(EXTRAS)/include

CONV_VERSION = 2.2.0

STD      = 98

TXT      = -DNO_SOURCE -DAUTO_DOC
DEFINES  =

COMPILE  = s.compile
FC       = $(COMPILE) -includes "$(MODDIR)" -defines "=$(DEFINES)" -O $(OPTIL) -optf "$(FFLAGS)" -bindir $(LCLPO) $(COMPF) $(FCOMPF) -src
FTNC     = $(COMPILE) -includes "$(MODDIR)" -defines "=$(DEFINES)" -P -src
FTNT     = $(CPPT) $(TXT)

CODEBETA =
#OBJ     = $(LIBDIR)/msg.o
#OBJ     = $(LIBDIR)/stubs.o
OBJ      =

MALIB    = cdf2ccc
LIBP     = $(LIBDIR)/lib$(MALIB).a

DDFUN90  = ddfun90
DIAGLIB  = diag_sq$(STD)_g
lNetCDF  = netcdff netcdf
RMNLIB   = rmn_015
UDUNITS  = udunits
VGDLIB   = descrip

LIBPATH  = $(LIBDIR) $(EXTRAS)/lib

# les programatheques externes


default:  $(LIBP)

.PRECIOUS: $(LIBP)

.DEFAULT:
	cvs checkout $@

.ftn.f:
	$(FTNC) $<

.ftn.o:
	$(FC) $<
	/bin/rm -f $*.f

.ftn.a:
	$(FC) $<
	$(AR) rv $@ $*.o
	/bin/rm -f $*.o

.f.o:
	$(FC) $<

.f.a:
	$(FC) $<
	$(AR) rv $@ $*.o
	/bin/rm -f $*.o

.o.a:
	$(AR) rv $@ $*.o

.ftn.txt:
	/bin/rm -f $(MANDIR)/$*.txt $*.txt
	$(FTNT) $< | sed s/'^[*]'/' '/ > $(MANDIR)/$*.txt
	rsync -av $(MANDIR)/$*.txt .
	if ( test "$(HOSTWEB):$(DIAGWEB)" != ":" ) then \
	rsync -av $(MANDIR)/$*.txt $(HOSTWEB):$(DIAGWEB) ; fi


FTNDECKS=\
	 decodate.ftn	encodate.ftn	enleve_bissextile2.ftn \
	 udparse.ftn	uduerr.ftn

FDECKS= \
	 affecte_attr.f	affecte_coord.f	affecte_dim.f	affecte_var.f \
	 attr_value.f	attribut_coord.f		attribut_var.f \
	 auxiliary_coordinate.f		 cal_date.f	clean_char.f	combline.f \
	 conv_unit.f	def_attr_projection.f		def_cccma.f \
	 def_date.f	def_dim_coord.f	def_level.f	def_name.f \
	 def_spectral_truncation.f	define_list.f	define_netcdf.f \
	 define_var.f	eval_lonlat.f	eval_xcyc.f \
	 get_attr.f	get_attribut.f	get_coord2.f	get_coordonne.f \
	 get_dim.f	get_name.f	get_string.f	get_topo.f \
	 get_var.f	handle_err.f \
	 hybrid_vertical_coordinate.f	init_dim.f	initialise.f	inq_file.f \
	 justifie_droite.f		les_arg.f	lire_arg.f	lire_cccvar.f \
	 minmaxchar.f	nblen.f	possibilite.f	put_attribut.f \
	 put_var.f	rdlatlon2.f		rdspec2.f	scanfile.f \
	 trier.f	valide_range.f	vers_cccma.f	vers_netcdf.f \
	 wrfile.f	wrlalo.f

EXTRAF= \
	cdf2ccc.f	decodate.f	encodate.f	enleve_bissextile2.f	udparse.f	uduerr.f


OBJECTS= \
	 affecte_attr.o	affecte_coord.o	affecte_dim.o	affecte_var.o \
	 attr_value.o	attribut_coord.o	attribut_var.o	auxiliary_coordinate.o \
	 cal_date.o			clean_char.o	combline.o \
	 conv_unit.o	decodate.o	def_attr_projection.o	def_cccma.o \
	 def_date.o	def_dim_coord.o	def_level.o	def_name.o \
	 def_spectral_truncation.o	define_list.o	define_netcdf.o	define_var.o \
	 encodate.o	enleve_bissextile2.o	eval_lonlat.o	eval_xcyc.o \
	 get_attr.o	get_attribut.o	get_coord2.o	get_coordonne.o \
	 get_dim.o	get_name.o	get_string.o	get_topo.o \
	 get_var.o	handle_err.o \
	 hybrid_vertical_coordinate.o	init_dim.o	initialise.o	inq_file.o \
	 justifie_droite.o	les_arg.o	lire_arg.o	lire_cccvar.o \
	 minmaxchar.o	nblen.o	possibilite.o	put_attribut.o \
	 put_var.o	rdlatlon2.o		rdspec2.o	scanfile.o \
	 trier.o		udparse.o	uduerr.o \
	 valide_range.o	vers_cccma.o	vers_netcdf.o	wrfile.o \
	 wrlalo.o

HDECKS= \
	$(INCLUDES)/cdf2ccc.h \
	$(INCLUDES)/dimmem.h \
	$(INCLUDES)/format.h \
	$(INCLUDES)/ibufmem.h \
	$(INCLUDES)/infomem.h \
	$(INCLUDES)/specmem.h \
	$(INCLUDES)/varmem.h \
	$(INCLUDES)/workmem.h \
	$(INCLUDES)/ztypmem.h

TXTDECKS =	cdf2ccc.txt



FICHIERS = $(FTNDECKS) $(FDECKS) $(HDECKS)

# Dependances des fichiers .txt et .o

cdf2ccc.txt:		cdf2ccc.ftn

$(LIBP)(affecte_attr.o):	affecte_attr.f $(INCLUDES)/cdf2ccc.h	$(INCLUDES)/workmem.h 
$(LIBP)(affecte_coord.o):	affecte_coord.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/infomem.h 
$(LIBP)(affecte_dim.o):		affecte_dim.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h 
$(LIBP)(affecte_var.o): 	affecte_var.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/infomem.h 
$(LIBP)(attr_value.o): 		attr_value.f	$(INCLUDES)/cdf2ccc.h
$(LIBP)(attribut_coord.o):	attribut_coord.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/format.h	$(INCLUDES)/infomem.h	$(INCLUDES)/varmem.h	$(INCLUDES)/workmem.h 
$(LIBP)(attribut_var.o): 	attribut_var.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/infomem.h 	$(INCLUDES)/varmem.h 
$(LIBP)(auxiliary_coordinate.o): 	auxiliary_coordinate.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/infomem.h	$(INCLUDES)/varmem.h	$(INCLUDES)/ztypmem.h
$(LIBP)(cal_date.o):			cal_date.f	$(INCLUDES)/cdf2ccc.h 
$(LIBP)(cdf2ccc.o):			cdf2ccc.ftn	$(INCLUDES)/cdf2ccc.h
$(LIBP)(combline.o): 		combline.f	$(INCLUDES)/cdf2ccc.h 
$(LIBP)(conv_unit.o): 		conv_unit.f	$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/infomem.h 
$(LIBP)(decodate.o):			decodate.ftn	$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/infomem.h  
$(LIBP)(def_attr_projection.o): def_attr_projection.f $(INCLUDES)/cdf2ccc.h $(INCLUDES)/ztypmem.h
$(LIBP)(def_cccma.o):		def_cccma.f	$(INCLUDES)/cdf2ccc.h 
$(LIBP)(def_dim_coord.o):	def_dim_coord.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/infomem.h	$(INCLUDES)/varmem.h		$(INCLUDES)/workmem.h 
$(LIBP)(def_level.o):		def_level.f	$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/dimmem.h 	$(INCLUDES)/infomem.h 	$(INCLUDES)/varmem.h		$(INCLUDES)/workmem.h 
$(LIBP)(def_spectral_truncation.o):	def_spectral_truncation.f	$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/infomem.h 
$(LIBP)(define_list.o):		define_list.f	$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/infomem.h 
$(LIBP)(define_netcdf.o):	define_netcdf.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/infomem.h 
$(LIBP)(define_var.o):		define_var.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/infomem.h 
$(LIBP)(encodate.o):			encodate.ftn	$(INCLUDES)/cdf2ccc.h 
$(LIBP)(enleve_bissextile2.o):	enleve_bissextile2.ftn	$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/format.h 
$(LIBP)(eval_lonlat.o):		eval_lonlat.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/infomem.h	$(INCLUDES)/varmem.h	$(INCLUDES)/ztypmem.h
$(LIBP)(eval_xcyc.o):		eval_xcyc.f	$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/infomem.h 	$(INCLUDES)/varmem.h 
$(LIBP)(get_attr.o):			get_attr.f		$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/workmem.h 
$(LIBP)(get_attribut.o):	get_attribut.f		$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/workmem.h 
$(LIBP)(get_coord2.o):		get_coord2.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/infomem.h	$(INCLUDES)/varmem.h 
$(LIBP)(get_coordonne.o):	get_coordonne.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/infomem.h	$(INCLUDES)/varmem.h 
$(LIBP)(get_dim.o):			get_dim.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h 
$(LIBP)(get_topo.o):			get_topo.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/varmem.h	$(INCLUDES)/workmem.h 
$(LIBP)(get_var.o):			get_var.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/workmem.h 
$(LIBP)(hybrid_vertical_coordinate.o):		hybrid_vertical_coordinate.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/infomem.h	$(INCLUDES)/varmem.h 
$(LIBP)(init_dim.o):			init_dim.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/infomem.h 
$(LIBP)(initialise.o):		initialise.f	$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/infomem.h
$(LIBP)(inq_file.o):			inq_file.f	$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/ztypmem.h
$(LIBP)(lire_arg.o):			lire_arg.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/ztypmem.h
$(LIBP)(lire_cccvar.o):		lire_cccvar.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/infomem.h	$(INCLUDES)/varmem.h	$(INCLUDES)/workmem.h 
$(LIBP)(possibilite.o):		possibilite.f	$(INCLUDES)/cdf2ccc.h 
$(LIBP)(put_attribut.o):	put_attribut.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/workmem.h 
$(LIBP)(put_var.o):			put_var.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/workmem.h 
$(LIBP)(rdlatlon2.o):		rdlatlon2.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/ibufmem.h	$(INCLUDES)/infomem.h	$(INCLUDES)/varmem.h	$(INCLUDES)/workmem.h	$(INCLUDES)/ztypmem.h 
$(LIBP)(rdspec2.o):			rdspec2.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h $(INCLUDES)/ibufmem.h	$(INCLUDES)/infomem.h $(INCLUDES)/infomem.h	$(INCLUDES)/specmem.h	$(INCLUDES)/varmem.h	$(INCLUDES)/workmem.h 
$(LIBP)(scanfile.o):			scanfile.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/format.h	$(INCLUDES)/varmem.h	$(INCLUDES)/workmem.h 
$(LIBP)(trier.o):				trier.f		$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/infomem.h	$(INCLUDES)/specmem.h $(INCLUDES)/ztypmem.h
$(LIBP)(valide_range.o):	valide_range.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/infomem.h 
$(LIBP)(vers_cccma.o):		vers_cccma.f	$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/dimmem.h	$(INCLUDES)/ibufmem.h	$(INCLUDES)/infomem.h	$(INCLUDES)/specmem.h	$(INCLUDES)/varmem.h	$(INCLUDES)/ztypmem.h
$(LIBP)(vers_netcdf.o):		vers_netcdf.f	$(INCLUDES)/cdf2ccc.h 	$(INCLUDES)/dimmem.h	$(INCLUDES)/ibufmem.h $(INCLUDES)/infomem.h	$(INCLUDES)/specmem.h	$(INCLUDES)/varmem.h	$(INCLUDES)/ztypmem.h
$(LIBP)(wrfile.o):			wrfile.f	$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/ibufmem.h	$(INCLUDES)/infomem.h	$(INCLUDES)/varmem.h	$(INCLUDES)/workmem.h 
$(LIBP)(wrlalo.o):			wrlalo.f		$(INCLUDES)/cdf2ccc.h	$(INCLUDES)/dimmem.h	$(INCLUDES)/ibufmem.h	$(INCLUDES)/infomem.h	$(INCLUDES)/varmem.h	$(INCLUDES)/workmem.h 


#LISTE DES CIBLES DU FICHIER recettes

obj: $(OBJECTS)
#Produire les fichiers objets (.o) pour tous les fichiers

sortirtout:
	extrait_tout

genlib: $(OBJECTS)
#Creer une programmatheque ayant pour nom $LIBP et incluant TOUS les fichiers objets
	$(AR) rv $(LIBP) $(OBJECTS)

libexp: sortirtout objexp
#Mettre tous les objets de l experience en cours dans la programmatheque $LIBP
	$(AR) rv $(LIBP) *.o

#Descriptions des recettes personnelles

initial:
# Setup des repertoires dependants de $(EC_ARCH)
	/bin/mkdir -p $(LIBDIR) $(MANDIR) $(MODDIR) $(BINDIR)
	s.locate --lib $(VGDLIB) 1> /dev/null || { echo -e \n PLS execute \". s.ssmuse.dot vgriddesc\" \n ; false ; }
	if [[ ! -f $(EXTRAS)/lib/libnetcdf.a ]]; then cd $(EXTRAS)/../../.. ; make all ; cd $(DIAGNOSTIQUE)/src ; fi

# Voir "Managing project with make" de O'Reilly & Associates, 1989,
# pages 51-53, pour les details de ce qui se passe dans ce qui suit
$(LIBP):	initial	$(FDECKS)	$(HDECKS)	\
	 $(LIBP)(affecte_attr.o)		\
	 $(LIBP)(affecte_coord.o)		\
	 $(LIBP)(affecte_dim.o)			\
	 $(LIBP)(affecte_var.o)			\
	 $(LIBP)(attr_value.o)			\
	 $(LIBP)(attribut_coord.o)		\
	 $(LIBP)(attribut_var.o)		\
	 $(LIBP)(auxiliary_coordinate.o)	\
	 $(LIBP)(cal_date.o)			\
	 $(LIBP)(clean_char.o)			\
	 $(LIBP)(combline.o) 			\
	 $(LIBP)(conv_unit.o)			\
	 $(LIBP)(decodate.o)			\
	 $(LIBP)(def_attr_projection.o)	\
	 $(LIBP)(def_cccma.o)			\
	 $(LIBP)(def_date.o)			\
	 $(LIBP)(def_dim_coord.o)		\
	 $(LIBP)(def_level.o)			\
	 $(LIBP)(def_name.o) 			\
	 $(LIBP)(def_spectral_truncation.o)	\
	 $(LIBP)(define_list.o)			\
	 $(LIBP)(define_netcdf.o)		\
	 $(LIBP)(define_var.o)			\
	 $(LIBP)(encodate.o)			\
	 $(LIBP)(enleve_bissextile2.o)		\
	 $(LIBP)(eval_lonlat.o)			\
	 $(LIBP)(eval_xcyc.o)			\
	 $(LIBP)(get_attr.o)			\
	 $(LIBP)(get_attribut.o)			\
	 $(LIBP)(get_coord2.o)			\
	 $(LIBP)(get_coordonne.o)		\
	 $(LIBP)(get_dim.o)			\
	 $(LIBP)(get_name.o)			\
	 $(LIBP)(get_string.o)			\
	 $(LIBP)(get_topo.o)			\
	 $(LIBP)(get_var.o)			\
	 $(LIBP)(handle_err.o)			\
	 $(LIBP)(hybrid_vertical_coordinate.o)	\
	 $(LIBP)(init_dim.o)			\
	 $(LIBP)(initialise.o)			\
	 $(LIBP)(inq_file.o)			\
	 $(LIBP)(justifie_droite.o)		\
	 $(LIBP)(les_arg.o)			\
	 $(LIBP)(lire_arg.o)			\
	 $(LIBP)(lire_cccvar.o)			\
	 $(LIBP)(minmaxchar.o)			\
	 $(LIBP)(nblen.o)			\
	 $(LIBP)(possibilite.o)			\
	 $(LIBP)(put_attribut.o)			\
	 $(LIBP)(put_var.o)			\
	 $(LIBP)(rdlatlon2.o)			\
	 $(LIBP)(rdspec2.o)			\
	 $(LIBP)(scanfile.o)			\
	 $(LIBP)(trier.o)			\
	 $(LIBP)(udparse.o)			\
	 $(LIBP)(uduerr.o)			\
	 $(LIBP)(valide_range.o)		\
	 $(LIBP)(vers_cccma.o)			\
	 $(LIBP)(vers_netcdf.o)			\
	 $(LIBP)(wrfile.o)			\
	 $(LIBP)(wrlalo.o)

document: cdf2ccc.ftn $(TXTDECKS)
# Creation de la documentation pour cdf2ccc

cdf2ccc: initial cdf2ccc.ftn
# Creation de l'executable cdf2ccc ailleurs que ailleurs que sous AIX
	$(AR) xv $(LIBP) lire_arg.o
	$(FC) $@.ftn -o $@ -libpath "$(LIBPATH)" \
	-libappl "$(MALIB) $(lNetCDF) $(UDUNITS) $(DIAGLIB) $(VGDLIB)" \
	-librmn $(RMNLIB) -obj "lire_arg.o $(OBJ)" $(CODEBETA) -libsys $(DDFUN90)
	/bin/rm lire_arg.o cdf2ccc.o ; /bin/mv $@ $(BINDIR)/cdf2ccc_$(CONV_VERSION) ;\
	cd $(BINDIR)  ; /bin/ln -f cdf2ccc_$(CONV_VERSION) cdf2rpn_$(CONV_VERSION)

cdf2rpn: cdf2ccc

AIX: initial cdf2ccc.ftn
# Creation de l'executable cdf2ccc pour architecture AIX
	$(AR) xv $(LIBP) lire_arg.o
	$(FC) cdf2ccc.ftn -o $@ -libpath "$(LIBPATH)" \
	-libappl "$(MALIB) $(lNetCDF) $(UDUNITS) $(DIAGLIB) $(VGDLIB)" \
	-librmn $(RMNLIB) -obj "lire_arg.o $(OBJ)" $(CODEBETA)
	/bin/rm lire_arg.o cdf2ccc.o ; /bin/mv $@ $(BINDIR)/cdf2ccc_$(CONV_VERSION) ;\
	cd $(BINDIR)  ; /bin/ln -f cdf2ccc_$(CONV_VERSION) cdf2rpn_$(CONV_VERSION)

clean:
# Faire un petit menage. On enleve tous les .o et les fichiers .f superflus
	/bin/rm -f *.o $(EXTRAF)
