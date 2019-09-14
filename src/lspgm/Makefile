.SUFFIXES :
.SUFFIXES : .ptn90 .ptn .txt .f90 .f .c .o .a

# Description: Makefile de LSPGM

# revision 4.162 2019/08/12 18:18  dugas
# - Utiliser libdescrip qui se trouve dans src/extras
#
# revision 4.161 2019/06/04 17:17  dugas
# - Passer a la version 6.4.5 du binaire r.diag (driver.Abs)
#
# Revision 4.160 2019/05/22 11:11  dugas
# - Ajouter la dependance a stats_signatures.mod pour AVERAGES
# - Ajouter la variable ARUFLAG (=U, par default)
#
# Revision 4.159 2019/05/15 15:15  dugas
# - Ajouter le programme AVERAGES (de MFV)
#
# Revision 4.158 2019/02/19 13:13  dugas
# - Enlever toutes re-definitions de MAKE qui pourraient inclure
#   une revision de la valeur de la variable DIAGNOSTIQUE. Cette
#   variable est deja correctement definie localement
#
# Revision 4.157 2019/02/12 14:14  dugas
# - Remplacer les instances de MAKE par des MYMAKE
#
# revision 4.156 2019/01/10 16:16  dugas
# - Specifier les dependances de la cible gshumh.o
# - Ajouter les dependances aux nouveaux comdecks
#   toc_gest_formats.cdk et set_a_b_formats.cdk
#
# revision 4.155 2018/11/16 17:17  dugas
# - Ajouter la dependance au module diag_toc a la cible windowf.o
#
# revision 4.154 2018/08/16 16:16  dugas
# - Passer a la version 6.4.4 du binaire r.diag (driver.Abs)
#
# revision 4.153 2018/07/21 13:13  dugas
# - Ajouter la dependance a zlabel.cdk pour la cible timavg.o
#
# revision 4.152 2018/01/28 17:17  dugas
# - Appel a command-change-date.sh lors de la generation du binaire r.diag
#
# revision 4.151 2017/11/02 15:15  dugas
# - Modifier/Completer les dependances aux modules de r.diag
#
# revision 4.150 2016/10/27 15:49  dugas
# - Le macro ENTETE passe de 28 a 32
# - Corriger la syntaxe de la liste de dependances des .o
# - Ajouter des dependances a toc_gest_1.cdk, toc_gest_2.cdk et set_a_b_from_AABB.cdk
#   pour DXDP, GPASG, GPDPD, GPHUM, GSAPL, GSAPZL, GSATHL, GSMSLPH, PGEN et VSINTH

ifeq ($(BASE_ARCH),)
$(error FATAL: BASE_ARCH is not defined or empty, ABORTING)
endif

ifeq "$(BASE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

SHELL   = /bin/bash

ARUFLAG = U
AR      = s.ar
ARR     = $(AR) rv$(ARUFLAG)

CPPT    = r.gppf -chop_bang -gpp -F
ZIP     = $(shell which zip)

DIAGNOSTIQUE = $(CURDIR)/../..

FFLAGS  =  

ENTETE  = 32
NP      = 1
OPTIL   = 1
STD     = 98

ARCH    = $(BASE_ARCH)

DIAGS   = $(DIAGNOSTIQUE)
BINDIR  = $(DIAGS)/bin/$(BASE_ARCH)
LIBDIR  = $(DIAGS)/lib/$(EC_ARCH)
MANDIR  = $(DIAGS)/man/pdoc
INCLUDE = $(DIAGS)/include
MODDIR  = $(INCLUDE)/$(EC_ARCH)

CODEBETA=
OBJ     =

BLAS    = blas
DDFUN90 = ddfun90
DIAGLIB = diag_sq$(STD)
GRAFLIB =
LAPACK  = lapack
PROGLIB = prog_sq$(STD)
RMNLIB  = rmn_016.3
VGDLIB  = descrip

LIBP    = $(LIBDIR)/lib$(PROGLIB).a
LIBS    = $(LIBDIR)/lib$(DIAGLIB).a
ZIPFILE = $(LIBDIR)/lspgm_STD$(STD).f.zip

DIAG_VERSION = 6.4.5

HOSTWEB =
DIAGWEB =

NOPLOT  = -DNOPLOT
DINTERP = -DNEW_INTERP
DEFINE  = $(NOPLOT) -Dnombre_de_taches=$(NP) -DVERY_NICE=oui -Dtaille_entete=$(ENTETE) $(DINTERP)
TXT     = -DNO_SOURCE -DAUTO_DOC -I $(INCLUDE)

COMPILE = s.compile
FC      = $(COMPILE) -includes "$(INCLUDE) $(MODDIR)" -defines "=$(DEFINE)" -O $(OPTIL) -optf="$(FFLAGS)" -src
CC      = $(COMPILE) -includes "$(INCLUDE)"           -defines "=$(DEFINE)" -O $(OPTIL) -src
FTNC    = $(COMPILE) -includes "$(INCLUDE)"           -defines "=$(DEFINE)" -P -src
FTNT    = $(CPPT) $(TXT)

default:  $(LIBP)

.PRECIOUS: $(LIBP)

.ptn.a:
	sed s/'      PROGRAM '/'      SUBROUTINE '/ $< > $*.ftn ;\
	$(FC) $*.ftn
	if [ -f $@ ]; then chmod 644 $@ ; fi
	$(ARR) $@ $*.o
	$(ZIP) $(ZIPFILE) $*.f
	/bin/rm -f $*.o $*.f $*.ftn

.ptn90.a:
	sed s/' PROGRAM $*'/' SUBROUTINE $*'/ $< > $*.ftn90 ;\
	$(FC) $*.ftn90
	if [ -f $@ ]; then chmod 644 $@ ; fi
	$(ARR) $@ $*.o
	$(ZIP) $(ZIPFILE) $*.f90
	/bin/rm -f $*.o $*.f90 $*.ftn90

.ptn.o:
	sed s/'      PROGRAM '/'      SUBROUTINE '/ $< > $*.ftn ;\
	$(FC) $*.ftn ;	/bin/rm -f $*.f $*.ftn

.ptn90.o:
	sed s/' PROGRAM $*'/' SUBROUTINE $*'/ $< > $*.ftn90 ;\
	$(FC) $*.ftn90 ;	/bin/rm -f $*.f90 $*.ftn90

.ptn.f:
	sed s/'      PROGRAM '/'      SUBROUTINE '/ $< > $*.ftn ;\
	$(FTNC) $*.ftn ; /bin/rm -f $*.ftn

.ptn90.f90:
	sed s/' PROGRAM $*'/' SUBROUTINE $*'/ $< > $*.ftn90 ;\
	$(FTNC) $*.ftn90 ; /bin/rm -f $*.ftn90

.ptn.txt:
	$(FTNT) $< | sed s/'^[*]'/' '/ > $(MANDIR)/$@ ;\
	$(DIAGNOSTIQUE)/bin/command-origin-document.sh $< $(MANDIR)/$@ ;\
	if ( test "$(HOSTWEB):$(DIAGWEB)" != ":" ) then \
	rsync -av $(MANDIR)/$@ $(HOSTWEB):$(DIAGWEB) ; fi

.ptn90.txt:
	$(FTNT) $< | sed s/'^[!]'/' '/ > $(MANDIR)/$@ ;\
	$(DIAGNOSTIQUE)/bin/command-origin-document.sh $< $(MANDIR)/$@ ;\
	if ( test "$(HOSTWEB):$(DIAGWEB)" != ":" ) then \
	rsync -av $(MANDIR)/$@ $(HOSTWEB):$(DIAGWEB) ; fi

.c.o:
	$(CC) $*.c


.f.a:

.f.o:

.f90.a:

.f90.o:


CDECKS =

HTMLDECKS = Change_Log.html	index.html	Tmpdir.html	Diag_Config.html

PTN90DECKS = \
	averages.ptn90		fstmdiag.ptn90

PTNDECKS = \
	addf.ptn	addtext.ptn	adjust.ptn	anomali.ptn	barx.ptn	bary.ptn \
	beta.ptn	betao.ptn	binach.ptn	binning.ptn	bins.ptn \
	calmap.ptn	calpres.ptn	cdiv.ptn	chabin.ptn	chklev.ptn \
	cmlt.ptn	cofagg.ptn	cofaps.ptn	complex.ptn	condef.ptn \
	conjmlt.ptn	constr.ptn	convert.ptn	corpar.ptn	cosf.ptn \
	critt.ptn	crvplot.ptn	ctemps.ptn	cwinds.ptn	delhat.ptn \
	delhato.ptn	delx.ptn	dely.ptn	dertim.ptn	diftim.ptn \
	divf.ptn	dxdp.ptn	enleve.ptn	eofrot.ptn	expone.ptn \
	extrac.ptn	fcats.ptn	fit.ptn 	flat.ptn	fmask.ptn	fpow.ptn \
	fptest.ptn	fsttri.ptn	ftest.ptn	fullcor.ptn	fullmix.ptn \
	gadd.ptn	gaussw.ptn	gdiv.ptn	ggacof.ptn	ggall.ptn \
	ggaps.ptn	ggatim.ptn	ggdev.ptn	ggdlat.ptn	ggdlon.ptn \
	gghems.ptn	gglob.ptn	ggplot.ptn	ggsign.ptn	ggstat.ptn \
	ggtrig.ptn	ggvecps.ptn	globavg.ptn	gmlt.ptn	gpasg.ptn \
	gpdpd.ptn	gphum.ptn	gpvtat.ptn	grset.ptn	gsapl.ptn \
	gsaptp.ptn	gsaptp.ptn	gsapzl.ptn	gsathl.ptn	gsathlt.ptn	gshumh.ptn \
	gsmslp.ptn	gsmslph.ptn	gsomgah.ptn	gsomega.ptn	gwtqd.ptn \
	hdplot.ptn	hovplot.ptn	ijlook.ptn	inttim.ptn	invntry.ptn \
	joiner.ptn	joinup.ptn	kef.ptn 	lagcor.ptn	laslon.ptn \
	legendre.ptn	levrev.ptn	linepc.ptn	linint.ptn	llagg.ptn \
	loge.ptn	lpprint.ptn	makezg.ptn	memcoef.ptn	mempred.ptn	mergef.ptn \
	meofvec.ptn	meriavg.ptn	minvers.ptn	miroir.ptn	mkmask.ptn	mksymm.ptn \
	mkunsy.ptn	mltf.ptn	mmult.ptn	modulus.ptn	tmmult.ptn \
	mslpr.ptn	mtrace.ptn	newnam.ptn	nspec.ptn 	parmsub.ptn \
	pcomp.ptn	pgen.ptn	plagalc.ptn	psdiv.ptn	psgrid.ptn \
	psmlt.ptn	randeof.ptn	rangrd.ptn	rcopy.ptn	relabl.ptn \
	repack.ptn	rgopr.ptn	rmax.ptn	rmlabl.ptn	rzonavg.ptn \
	select.ptn	select2.ptn	shiftf.ptn	signag.ptn	sinf.ptn \
	slice.ptn	smooth.ptn	smspec.ptn	solflux.ptn	sortf.ptn	spamp.ptn \
	spblow.ptn	spcalpa.ptn	spcdns.ptn	spconj.ptn	spcsum.ptn \
	spdlat.ptn	spdlon.ptn	spdlsq.ptn	spectr.ptn	spflux.ptn \
	spgdev.ptn	splapl.ptn	spldiss.ptn	splinv.ptn	splook.ptn \
	spnegm.ptn	spplot.ptn	spreim.ptn	spresol.ptn	spslope.ptn \
	spsmoo.ptn	spvdgx.ptn	spzdev.ptn	spzxmlt.ptn	sqroot.ptn \
	square.ptn	subarea.ptn	subf.ptn	t2test.ptn	tfilter.ptn \
	timagg.ptn	timavg.ptn	timcov.ptn	timdev.ptn	timegen.ptn \
	transpos.ptn	tsafc.ptn	tscheck.ptn	tstep.ptn	ttesta.ptn \
	tteste.ptn	txtplot.ptn	unbinn.ptn	unmask.ptn	unseas.ptn \
	varspec.ptn	vbasd.ptn	veint.ptn	vertinc.ptn	vminmax.ptn	vpint.ptn \
	vpintl.ptn	vpintv.ptn	vsint.ptn	vsinth.ptn	w2stat.ptn \
	windowf.ptn	xappend.ptn	xfind.ptn	xlin.ptn	xmplot.ptn \
	xpand.ptn	xplot.ptn	xsave.ptn	xtrac.ptn	xylin.ptn \
	ytox.ptn	zadlev.ptn	zdhdiv.ptn	zonavg.ptn	zoncov.ptn \
	zondev.ptn	zonmod.ptn	zxaglb.ptn	zxdev.ptn	zxdevp.ptn \
	zxdlat.ptn	zxdlatf.ptn	zxgam.ptn	zxint.ptn	zxplot.ptn \
	zxprint.ptn	zxpsi.ptn	zxtrig.ptn	zxwave.ptn	zxmav.ptn

TXTDECKS = \
	addf.txt	addtext.txt	adjust.txt	anomali.txt	averages.txt	barx.txt \
	bary.txt	beta.txt	betao.txt	binach.txt	binning.txt	bins.txt \
	calmap.txt	calpres.txt	cdiv.txt	chabin.txt	chklev.txt \
	cmlt.txt	cofagg.txt	cofaps.txt	complex.txt	condef.txt \
	conjmlt.txt	constr.txt	convert.txt	corpar.txt	cosf.txt \
	critt.txt	crvplot.txt	ctemps.txt	cwinds.txt	delhat.txt \
	delhato.txt	delx.txt	dely.txt	dertim.txt	diftim.txt \
	divf.txt	dxdp.txt	enleve.txt	eofrot.txt	expone.txt \
	extrac.txt	fcats.txt	fit.txt 	flat.txt	fmask.txt	fpow.txt \
	fptest.txt	fstmdiag.txt	fsttri.txt	ftest.txt	fullcor.txt	fullmix.txt \
	gadd.txt	gaussw.txt	gdiv.txt	ggacof.txt	ggall.txt \
	ggaps.txt	ggatim.txt	ggdev.txt	ggdlat.txt	ggdlon.txt \
	gghems.txt	gglob.txt	ggplot.txt	ggsign.txt	ggstat.txt \
	ggtrig.txt	ggvecps.txt	globavg.txt	gmlt.txt	gpasg.txt \
	gpdpd.txt	gphum.txt	gpvtat.txt	grset.txt	gsapl.txt \
	gsaptp.txt	gsaptp.txt	gsapzl.txt	gsathl.txt	gsathlt.txt	gshumh.txt \
	gsmslp.txt	gsmslph.txt	gsomgah.txt	gsomega.txt	gwtqd.txt \
	hdplot.txt	hovplot.txt	ijlook.txt	inttim.txt	invntry.txt \
	joiner.txt	joinup.txt	kef.txt 	lagcor.txt	laslon.txt \
	legendre.txt	levrev.txt	linepc.txt	linint.txt	llagg.txt \
	loge.txt	lpprint.txt	makezg.txt	memcoef.txt	mempred.txt	mergef.txt \
	meofvec.txt	meriavg.txt	minvers.txt	miroir.txt	mkmask.txt	mksymm.txt \
	mkunsy.txt	mltf.txt	mmult.txt	modulus.txt	tmmult.txt \
	mslpr.txt	mtrace.txt	newnam.txt	nspec.txt 	parmsub.txt \
	pcomp.txt	pgen.txt	plagalc.txt	psdiv.txt	psgrid.txt \
	psmlt.txt	randeof.txt	rangrd.txt	rcopy.txt	relabl.txt \
	repack.txt	rgopr.txt	rmax.txt	rmlabl.txt	rzonavg.txt \
	select.txt	select2.txt	shiftf.txt	signag.txt	sinf.txt \
	slice.txt	smooth.txt	smspec.txt	solflux.txt	sortf.txt	spamp.txt \
	spblow.txt	spcalpa.txt	spcdns.txt	spconj.txt	spcsum.txt \
	spdlat.txt	spdlon.txt	spdlsq.txt	spectr.txt	spflux.txt \
	spgdev.txt	splapl.txt	spldiss.txt	splinv.txt	splook.txt \
	spnegm.txt	spplot.txt	spreim.txt	spresol.txt	spslope.txt \
	spsmoo.txt	spvdgx.txt	spzdev.txt	spzxmlt.txt	sqroot.txt \
	square.txt	subarea.txt	subf.txt	t2test.txt	tfilter.txt \
	timagg.txt	timavg.txt	timcov.txt	timdev.txt	timegen.txt \
	transpos.txt	tsafc.txt	tscheck.txt	tstep.txt	ttesta.txt \
	tteste.txt	txtplot.txt	unbinn.txt	unmask.txt	unseas.txt \
	varspec.txt	vbasd.txt	veint.txt	vertinc.txt	vminmax.txt	vpint.txt \
	vpintl.txt	vpintv.txt	vsint.txt	vsinth.txt	w2stat.txt \
	windowf.txt	xappend.txt	xfind.txt	xlin.txt	xmplot.txt \
	xpand.txt	xplot.txt	xsave.txt	xtrac.txt	xylin.txt \
	ytox.txt	zadlev.txt	zdhdiv.txt	zonavg.txt	zoncov.txt \
	zondev.txt	zonmod.txt	zxaglb.txt	zxdev.txt	zxdevp.txt \
	zxdlat.txt	zxdlatf.txt	zxgam.txt	zxint.txt	zxplot.txt \
	zxprint.txt	zxpsi.txt	zxtrig.txt	zxwave.txt	zxmav.txt

#LISTE DES MODULES AYANT DES DEPENDANCES

$(LIBP)(addf.o)    : addf.ptn \
	$(INCLUDE)/diagmacros.cdk
$(LIBP)(adjust.o)  : adjust.ptn \
	$(INCLUDE)/machtype.cdk
$(LIBP)(averages.o) : averages.ptn90 \
	$(MODDIR)/averages_common.mod \
	$(MODDIR)/stats_signatures.mod
$(LIBP)(binach.o)  : binach.ptn \
	$(INCLUDE)/machtype.cdk
$(LIBP)(chabin.o)  : chabin.ptn \
	$(INCLUDE)/machtype.cdk
$(LIBP)(convert.o) : convert.ptn \
	$(INCLUDE)/machtype.cdk
$(LIBP)(cosf.o)    : cosf.ptn \
	$(INCLUDE)/diagmacros.cdk
$(LIBP)(crvplot.o) : crvplot.ptn \
	$(INCLUDE)/char2int.cdk
$(LIBP)(dertim.o)  : dertim.ptn \
	$(INCLUDE)/diagmacros.cdk  $(INCLUDE)/machtype.cdk 
$(LIBP)(diftim.o)  : diftim.ptn \
	$(INCLUDE)/zlabel.cdk
$(LIBP)(dxdp.o)    : dxdp.ptn \
	$(INCLUDE)/toc_gest_1.cdk	$(INCLUDE)/toc_gest_2.cdk \
	$(INCLUDE)/toc_gest_formats.cdk	$(INCLUDE)/set_a_b_from_AABB.cdk \
	$(INCLUDE)/set_a_b_formats.cdk	$(MODDIR)/diag_toc.mod
$(LIBP)(expone.o)  : expone.ptn \
	$(INCLUDE)/diagmacros.cdk
$(LIBP)(fcats.o)   : fcats.ptn \
	$(INCLUDE)/machtype.cdk
$(LIBP)(fstmdiag.o) : fstmdiag.ptn90 \
	$(MODDIR)/diag_convert_ip123.mod \
	$(MODDIR)/stats_signatures.mod
$(LIBP)(fsttri.o) : fsttri.ptn \
	$(MODDIR)/diag_convert_ip123.mod
$(LIBP)(gadd.o)    : gadd.ptn \
	$(INCLUDE)/diagmacros.cdk
$(LIBP)(ggatim.o)  : ggatim.ptn \
	$(INCLUDE)/diagmacros.cdk
$(LIBP)(ggplot.o)  : ggplot.ptn \
	$(INCLUDE)/char2int.cdk
$(LIBP)(gpasg.o) : gpasg.ptn \
	$(INCLUDE)/toc_gest_1.cdk	$(INCLUDE)/toc_gest_formats.cdk \
	$(INCLUDE)/set_a_b_from_AABB.cdk	$(INCLUDE)/set_a_b_formats.cdk \
	$(MODDIR)/diag_toc.mod
$(LIBP)(gpdpd.o)   : gpdpd.ptn \
	$(INCLUDE)/consdia.cdk	$(INCLUDE)/dintern.cdk \
	$(INCLUDE)/fintern.cdk	$(INCLUDE)/set_a_b_from_AABB.cdk \
	$(INCLUDE)/set_a_b_formats.cdk	$(INCLUDE)/toc_gest_1.cdk \
	$(INCLUDE)/toc_gest_2.cdk	$(INCLUDE)/toc_gest_formats.cdk \
	$(MODDIR)/diag_toc.mod
$(LIBP)(gphum.o)   : gphum.ptn \
	$(INCLUDE)/consdia.cdk $(INCLUDE)/dintern.cdk \
	$(INCLUDE)/fintern.cdk	$(INCLUDE)/toc_gest_1.cdk \
	$(INCLUDE)/toc_gest_2.cdk	$(INCLUDE)/set_a_b_from_AABB.cdk \
	$(INCLUDE)/set_a_b_formats.cdk	$(MODDIR)/diag_toc.mod
$(LIBP)(gpvtat.o)  : gpvtat.ptn \
	$(INCLUDE)/consdia.cdk $(INCLUDE)/dintern.cdk \
	$(INCLUDE)/fintern.cdk
$(LIBP)(gsapl.o) : gsapl.ptn \
	$(INCLUDE)/toc_gest_1.cdk	$(INCLUDE)/toc_gest_2.cdk \
	$(INCLUDE)/toc_gest_formats.cdk	$(INCLUDE)/set_a_b_from_AABB.cdk \
	$(INCLUDE)/set_a_b_formats.cdk	$(MODDIR)/diag_toc.mod
$(LIBP)(gsapzl.o) : gsapzl.ptn \
	$(INCLUDE)/toc_gest_1.cdk	$(INCLUDE)/toc_gest_2.cdk \
	$(INCLUDE)/toc_gest_formats.cdk	$(INCLUDE)/set_a_b_from_AABB.cdk \
	$(INCLUDE)/set_a_b_formats.cdk	$(MODDIR)/diag_toc.mod
$(LIBP)(gsathl.o) : gsathl.ptn \
	$(INCLUDE)/toc_gest_1.cdk	$(INCLUDE)/toc_gest_2.cdk \
	$(INCLUDE)/toc_gest_formats.cdk	$(INCLUDE)/set_a_b_from_AABB.cdk \
	$(INCLUDE)/set_a_b_formats.cdk	$(MODDIR)/diag_toc.mod
$(LIBP)(gshumh.o) : gshumh.ptn \
	$(INCLUDE)/consdia.cdk $(INCLUDE)/dintern.cdk \
	$(INCLUDE)/fintern.cdk	$(INCLUDE)/toc_gest_1.cdk \
	$(INCLUDE)/toc_gest_formats.cdk	$(MODDIR)/diag_toc.mod
$(LIBP)(gsmslph.o) : gsmslph.ptn \
	$(INCLUDE)/toc_gest_1.cdk	$(INCLUDE)/toc_gest_2.cdk \
	$(INCLUDE)/toc_gest_formats.cdk	$(INCLUDE)/set_a_b_from_AABB.cdk \
	$(INCLUDE)/set_a_b_formats.cdk	$(MODDIR)/diag_toc.mod
$(LIBP)(hovplot.o) : hovplot.ptn \
	$(INCLUDE)/char2int.cdk
$(LIBP)(inttim.o)  : inttim.ptn \
	$(INCLUDE)/diagmacros.cdk $(INCLUDE)/machtype.cdk 
$(LIBP)(invntry.o) : invntry.ptn \
	$(INCLUDE)/machtype.cdk
$(LIBP)(joiner.o)  : joiner.ptn \
	$(INCLUDE)/machtype.cdk
$(LIBP)(meofvec.o) : meofvec.ptn \
	$(INCLUDE)/diagmacros.cdk $(INCLUDE)/machtype.cdk 
$(LIBP)(minvers.o) : minvers.ptn \
	$(INCLUDE)/diagmacros.cdk $(INCLUDE)/machtype.cdk 
$(LIBP)(mmult.o)   : mmult.ptn \
	$(INCLUDE)/diagmacros.cdk
$(LIBP)(mtrace.o)  : mtrace.ptn \
	$(INCLUDE)/diagmacros.cdk
$(LIBP)(pgen.o) : pgen.ptn \
	$(INCLUDE)/toc_gest_1.cdk	$(INCLUDE)/toc_gest_2.cdk \
	$(INCLUDE)/toc_gest_formats.cdk	$(INCLUDE)/set_a_b_from_AABB.cdk \
	$(INCLUDE)/set_a_b_formats.cdk	$(MODDIR)/diag_toc.mod
$(LIBP)(psdiv.o)   : psdiv.ptn \
	$(INCLUDE)/diagmacros.cdk
$(LIBP)(psmlt.o)   : psmlt.ptn \
	$(INCLUDE)/diagmacros.cdk
$(LIBP)(select.o)  : select.ptn \
	$(INCLUDE)/zlabel.cdk \
	$(MODDIR)/diag_toc.mod
$(LIBP)(select2.o) : select2.ptn \
	$(INCLUDE)/zlabel.cdk \
	$(MODDIR)/diag_toc.mod
$(LIBP)(shiftf.o)  : shiftf.ptn \
	$(INCLUDE)/diagmacros.cdk $(INCLUDE)/machtype.cdk 
$(LIBP)(sinf.o)    : sinf.ptn \
	$(INCLUDE)/diagmacros.cdk
$(LIBP)(sortf.o)   : sortf.ptn \
	$(INCLUDE)/diagmacros.cdk $(INCLUDE)/machtype.cdk 
$(LIBP)(spplot.o) : spplot.ptn \
	$(INCLUDE)/char2int.cdk
$(LIBP)(spvdgx.o)  : spvdgx.ptn \
	$(INCLUDE)/machtype.cdk 
$(LIBP)(tfilter.o) : tfilter.ptn \
	$(INCLUDE)/diagmacros.cdk $(INCLUDE)/machtype.cdk 
$(LIBP)(timavg.o)  : timavg.ptn \
	$(INCLUDE)/diagmacros.cdk \
	$(INCLUDE)/machtype.cdk \
	$(INCLUDE)/zlabel.cdk \
	$(MODDIR)/diag_convert_ip123.mod \
	$(MODDIR)/stats_signatures.mod
$(LIBP)(timcov.o)  : timcov.ptn \
	$(INCLUDE)/diagmacros.cdk $(INCLUDE)/machtype.cdk \
	$(MODDIR)/stats_signatures.mod
$(LIBP)(timdev.o)  : timdev.ptn \
	$(INCLUDE)/diagmacros.cdk $(INCLUDE)/machtype.cdk 
$(LIBP)(timegen.o) : timegen.ptn \
	$(INCLUDE)/diagmacros.cdk
$(LIBP)(tsafc.o)   : tsafc.ptn \
	$(INCLUDE)/machtype.cdk
$(LIBP)(tscheck.o) : tscheck.ptn \
	$(INCLUDE)/machtype.cdk
$(LIBP)(tstep.o)   : tstep.ptn \
	$(INCLUDE)/diagmacros.cdk $(INCLUDE)/machtype.cdk 
$(LIBP)(ttesta.o)  : ttesta.ptn \
	$(INCLUDE)/zlabel.cdk
$(LIBP)(tteste.o)  : tteste.ptn \
	$(INCLUDE)/zlabel.cdk
$(LIBP)(txtplot.o) : txtplot.ptn \
	$(INCLUDE)/char2int.cdk
$(LIBP)(vsinth.o) : vsinth.ptn \
	$(INCLUDE)/toc_gest_1.cdk	$(INCLUDE)/toc_gest_2.cdk \
	$(INCLUDE)/toc_gest_formats.cdk	$(INCLUDE)/set_a_b_from_AABB.cdk \
	$(INCLUDE)/set_a_b_formats.cdk	$(MODDIR)/diag_toc.mod
$(LIBP)(windowf.o) : windowf.ptn	$(MODDIR)/diag_toc.mod
$(LIBP)(xappend.o) : xappend.ptn \
	$(INCLUDE)/machtype.cdk
$(LIBP)(xfind.o)   : xfind.ptn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk
$(LIBP)(xsave.o)   : xsave.ptn \
	$(INCLUDE)/machtype.cdk $(INCLUDE)/zlabel.cdk
$(LIBP)(zoncov.o)  : zoncov.ptn \
	$(INCLUDE)/diagmacros.cdk
$(LIBP)(zxplot.o) : zxplot.ptn \
	$(INCLUDE)/char2int.cdk

#LISTE DES CIBLES DU FICHIER recettes

initial:
# Setup des repertoires dependants de $(EC_ARCH) et de $(VGDLIB)
	/bin/mkdir -p $(LIBDIR) $(MANDIR) $(MODDIR) $(BINDIR)
#	s.locate --lib $(VGDLIB) 1> /dev/null || { echo -e \n PLS execute \". s.ssmuse.dot vgriddesc\" \n ; false ; }
	if [[ ! -f $(LIBDIR)/libdescrip.a || -z "$(VGDLIB)" ]]; then \
		cd $(DIAGNOSTIQUE)/src/extras/descrip ; $(MAKE) ARUFLAG=$(ARUFLAG) ; fi

# Voir "Managing project with make" de O'Reilly & Associates, 1989,
# pages 51-53, pour les details de ce qui se passe dans ce qui suit
$(LIBP):	initial	$(PTNDECKS)	$(PTN90DECKS) \
	$(LIBP)(addf.o) 	$(LIBP)(addtext.o)	$(LIBP)(adjust.o)	$(LIBP)(anomali.o)	$(LIBP)(averages.o) \
	$(LIBP)(barx.o) 	$(LIBP)(bary.o)	 	$(LIBP)(beta.o)	 	$(LIBP)(betao.o)	$(LIBP)(binach.o) \
	$(LIBP)(binning.o)	$(LIBP)(bins.o) 	$(LIBP)(calmap.o)		$(LIBP)(calpres.o) \
	$(LIBP)(cdiv.o) 	$(LIBP)(chabin.o)		$(LIBP)(chklev.o)		$(LIBP)(cmlt.o) \
	$(LIBP)(cofagg.o)	$(LIBP)(cofaps.o)		$(LIBP)(complex.o)	$(LIBP)(condef.o) \
	$(LIBP)(conjmlt.o)	$(LIBP)(convert.o)	$(LIBP)(corpar.o)	$(LIBP)(cosf.o) \
	$(LIBP)(critt.o)	$(LIBP)(crvplot.o)	$(LIBP)(ctemps.o)	$(LIBP)(cwinds.o) \
	$(LIBP)(delhat.o)	$(LIBP)(delhato.o)	$(LIBP)(dertim.o)	$(LIBP)(delx.o) \
	$(LIBP)(dely.o) 	$(LIBP)(diftim.o)	$(LIBP)(divf.o) 	$(LIBP)(dxdp.o) \
	$(LIBP)(enleve.o)	$(LIBP)(eofrot.o)	$(LIBP)(expone.o)	$(LIBP)(extrac.o) \
	$(LIBP)(fcats.o)	$(LIBP)(fit.o)  	$(LIBP)(flat.o) 	$(LIBP)(fmask.o)	$(LIBP)(fpow.o) \
	$(LIBP)(fptest.o)	$(LIBP)(fstmdiag.o)	$(LIBP)(fsttri.o)	$(LIBP)(ftest.o)	$(LIBP)(fullcor.o) \
	$(LIBP)(fullmix.o)	$(LIBP)(gadd.o) 	$(LIBP)(gaussw.o)	$(LIBP)(gdiv.o) \
	$(LIBP)(gmlt.o) 	$(LIBP)(ggacof.o)	$(LIBP)(ggall.o)	$(LIBP)(ggaps.o) \
	$(LIBP)(ggatim.o)	$(LIBP)(ggdev.o)	$(LIBP)(ggdlat.o)	$(LIBP)(ggdlon.o) \
	$(LIBP)(gghems.o)	$(LIBP)(gglob.o)	$(LIBP)(ggplot.o)	$(LIBP)(ggsign.o) \
	$(LIBP)(ggstat.o)	$(LIBP)(ggtrig.o)	$(LIBP)(ggvecps.o)	$(LIBP)(globavg.o) \
	$(LIBP)(gpasg.o)	$(LIBP)(gpdpd.o)	$(LIBP)(gphum.o)	$(LIBP)(gpvtat.o) \
	$(LIBP)(grset.o)	$(LIBP)(gsapl.o)	$(LIBP)(gsaptp.o)	$(LIBP)(gsapzl.o) \
	$(LIBP)(gsathl.o)	$(LIBP)(gsathlt.o)	$(LIBP)(gshumh.o)	$(LIBP)(gsmslp.o) \
	$(LIBP)(gsmslph.o)	$(LIBP)(gsomgah.o)	$(LIBP)(gsomega.o)	$(LIBP)(gwtqd.o) \
	$(LIBP)(hdplot.o)	$(LIBP)(hovplot.o)	$(LIBP)(ijlook.o)	$(LIBP)(inttim.o) \
	$(LIBP)(invntry.o)	$(LIBP)(joiner.o)	$(LIBP)(joinup.o)	$(LIBP)(kef.o) \
	$(LIBP)(lagcor.o)	$(LIBP)(laslon.o)	$(LIBP)(legendre.o)	$(LIBP)(levrev.o) \
	$(LIBP)(linepc.o)	$(LIBP)(linint.o)	$(LIBP)(llagg.o)	$(LIBP)(loge.o) \
	$(LIBP)(lpprint.o)	$(LIBP)(makezg.o)	$(LIBP)(memcoef.o)	$(LIBP)(mempred.o)	$(LIBP)(mergef.o) \
	$(LIBP)(meofvec.o)	$(LIBP)(meriavg.o)	$(LIBP)(minvers.o)	$(LIBP)(miroir.o)	$(LIBP)(mkmask.o) \
	$(LIBP)(mksymm.o)	$(LIBP)(mkunsy.o)	$(LIBP)(mltf.o) 	$(LIBP)(mmult.o) \
	$(LIBP)(modulus.o)	$(LIBP)(tmmult.o)	$(LIBP)(mslpr.o)	$(LIBP)(mtrace.o) \
	$(LIBP)(newnam.o)	$(LIBP)(nspec.o)	$(LIBP)(parmsub.o)	$(LIBP)(pcomp.o) \
	$(LIBP)(pgen.o) 	$(LIBP)(plagalc.o)	$(LIBP)(psdiv.o)	$(LIBP)(psgrid.o) \
	$(LIBP)(psmlt.o)	$(LIBP)(randeof.o)	$(LIBP)(rangrd.o)	$(LIBP)(rcopy.o) \
	$(LIBP)(relabl.o)	$(LIBP)(repack.o)	$(LIBP)(rgopr.o)	$(LIBP)(rmax.o) \
	$(LIBP)(rmlabl.o)	$(LIBP)(rzonavg.o)	$(LIBP)(select.o)	$(LIBP)(select2.o) \
	$(LIBP)(shiftf.o)	$(LIBP)(signag.o)	$(LIBP)(sinf.o) 	$(LIBP)(slice.o) \
	$(LIBP)(smooth.o)	$(LIBP)(smspec.o)	$(LIBP)(solflux.o)	$(LIBP)(sortf.o)	$(LIBP)(spamp.o) \
	$(LIBP)(spblow.o)	$(LIBP)(spcalpa.o)	$(LIBP)(spcdns.o)	$(LIBP)(spconj.o) \
	$(LIBP)(spcsum.o)	$(LIBP)(spdlat.o)	$(LIBP)(spdlon.o)	$(LIBP)(spdlsq.o) \
	$(LIBP)(spectr.o)	$(LIBP)(spflux.o)	$(LIBP)(spgdev.o)	$(LIBP)(splapl.o) \
	$(LIBP)(spldiss.o)	$(LIBP)(splinv.o)	$(LIBP)(splook.o)	$(LIBP)(spnegm.o) \
	$(LIBP)(spplot.o)	$(LIBP)(spreim.o)	$(LIBP)(spresol.o)	$(LIBP)(spslope.o) \
	$(LIBP)(spsmoo.o)	$(LIBP)(spvdgx.o)	$(LIBP)(spzdev.o)	$(LIBP)(spzxmlt.o) \
	$(LIBP)(sqroot.o)	$(LIBP)(square.o)	$(LIBP)(subarea.o)	$(LIBP)(subf.o) \
	$(LIBP)(t2test.o)	$(LIBP)(tfilter.o)	$(LIBP)(timagg.o)	$(LIBP)(timavg.o) \
	$(LIBP)(timcov.o)	$(LIBP)(timdev.o)	$(LIBP)(timegen.o)	$(LIBP)(transpos.o) \
	$(LIBP)(tsafc.o)	$(LIBP)(tscheck.o)	$(LIBP)(tstep.o)	$(LIBP)(ttesta.o) \
	$(LIBP)(tteste.o)	$(LIBP)(txtplot.o)	$(LIBP)(unbinn.o)	$(LIBP)(unmask.o) \
	$(LIBP)(unseas.o)	$(LIBP)(varspec.o)	$(LIBP)(vbasd.o)	$(LIBP)(veint.o) \
	$(LIBP)(vertinc.o)	$(LIBP)(vminmax.o)	$(LIBP)(vpint.o)	$(LIBP)(vpintl.o)	$(LIBP)(vpintv.o) \
	$(LIBP)(vsint.o)	$(LIBP)(vsinth.o)	$(LIBP)(w2stat.o)	$(LIBP)(windowf.o) \
	$(LIBP)(xappend.o)	$(LIBP)(xfind.o)	$(LIBP)(xlin.o) 	$(LIBP)(xmplot.o) \
	$(LIBP)(xpand.o)	$(LIBP)(xplot.o)	$(LIBP)(xsave.o)	$(LIBP)(xtrac.o) \
	$(LIBP)(xylin.o)	$(LIBP)(ytox.o)		$(LIBP)(zadlev.o)	$(LIBP)(zdhdiv.o) \
	$(LIBP)(zonavg.o)	$(LIBP)(zoncov.o)	$(LIBP)(zondev.o)	$(LIBP)(zonmod.o) \
	$(LIBP)(zxaglb.o)	$(LIBP)(zxdev.o)	$(LIBP)(zxdevp.o)	$(LIBP)(zxdlat.o) \
	$(LIBP)(zxdlatf.o)	$(LIBP)(zxgam.o)	$(LIBP)(zxint.o)	$(LIBP)(zxmav.o) \
	$(LIBP)(zxplot.o)	$(LIBP)(zxprint.o)	$(LIBP)(zxpsi.o)	$(LIBP)(zxtrig.o) \
	$(LIBP)(zxwave.o)

document: driver.ptn info.lspgm $(TXTDECKS)
#Section de creation de la documentation pour lspgm
	if [[ "$(HOSTWEB):$(DIAGWEB)" != ":" ]] ; then \
	rsync -av $(HTMLDECKS) $(HOSTWEB):$(DIAGWEB) ; fi
	$(MAKE) driver.txt

driver: driver.ptn $(LIBP) $(LIBS)
#Section de creation du binaire proprement dit
	$(CPPT) $(DEFINE) $@.ptn > $@.f
	$(DIAGNOSTIQUE)/bin/command-change-date.sh $@.f ; \
	$(FC) $@.f  -o $@ -libpath $(DIAGS)/lib/$(EC_ARCH) \
	-libappl "$(PROGLIB) $(DIAGLIB) $(VGDLIB) $(GRAFLIB)" \
	-librmn $(RMNLIB) $(CODEBETA) -obj "$(OBJ)" \
	-libsys "$(DDFUN90) $(LAPACK) $(BLAS)"
	if [ -x $@ ] ; then \
		/bin/mv driver $(BINDIR)/r.diag_$(DIAG_VERSION)_$(COMP_ARCH) ; \
		ssh localhost ldd $(BINDIR)/r.diag_$(DIAG_VERSION)_$(COMP_ARCH) | \
		grep -i "not found" && exit 1 ; /bin/rm $@.f $@.o ; \
		touch -t 199001010000 driver && chmod 755 driver ; \
	fi

info.lspgm: driver.ptn
	for i in a b c d e f g; do \
	/bin/rm -f $(MANDIR)/lspgm_$$i.txt $(MANDIR)/lspgm_$$i.html; \
	$(CPPT) -DNO_SOURCE -DMORE_INFO driver.ptn | grep "^[*]$$i" | sed "s/^[*]$$i/ /" > lspgm_$$i.txt; \
	/bin/cat lspgm_$$i.txt | info2html > lspgm_$$i.html ; \
	rsync -av lspgm_$$i.txt lspgm_$$i.html $(MANDIR) ; \
	if [[ "$(HOSTWEB):$(DIAGWEB)" != ":" ]] ; then \
	rsync -av lspgm_$$i.txt lspgm_$$i.html $(HOSTWEB):$(DIAGWEB) ; fi \
	done

AIX: #Section appropriee pour les fichiers STD saveur 1998 sur IBM AIX
	$(MAKE) driver OBJ="$(OBJ)" RMNLIB=$(RMNLIB)  DDFUN90=" " \
	NOPLOT="-DNOPLOT" GRAFLIB=" " VGDLIB=$(VGDLIB) \
	MANDIR=$(MANDIR)

AIX-powerpc7: #Section appropriee pour les fichiers STD saveur 1998 sur AIX-powerpc7
	$(MAKE) driver OBJ="$(OBJ)" RMNLIB=$(RMNLIB) DDFUN90=" " \
	NOPLOT=$(NOPLOT) GRAFLIB=$(GRAFLIB) LAPACK=lapack_340 \
	VGDLIB=$(VGDLIB) MANDIR=$(MANDIR)

Linux:
#Section appropriee pour les fichiers STD saveur 1998 sur Linux INTEL
	$(MAKE) driver OBJ="$(OBJ)" RMNLIB=$(RMNLIB) MANDIR=$(MANDIR) \
	NOPLOT=$(NOPLOT) GRAFLIB=$(GRAFLIB) VGDLIB=$(VGDLIB)

Linux_x86-64:
#Section appropriee pour les fichiers STD saveur 1998 sur Linux INTEL 64 bits
	$(MAKE) driver OBJ="$(OBJ)" RMNLIB=$(RMNLIB) MANDIR=$(MANDIR) \
	NOPLOT=$(NOPLOT) VGDLIB=$(VGDLIB) GRAFLIB=$(GRAFLIB)

NoPlot:
#Section appropriee pour les fichiers STD saveur 1998 sur Linux 64 bits
# mais sans les modules graphiques
	$(MAKE) driver OBJ="$(OBJ)" RMNLIB=$(RMNLIB) MANDIR=$(MANDIR) \
	NOPLOT="-DNOPLOT" GRAFLIB=" " VGDLIB=$(VGDLIB)

clean:
# Faire un petit menage. On enleve tous les .o
	/bin/rm -f *.o
