.SUFFIXES :
.SUFFIXES : .F90 .ptn90 .ptn .ftn .txt

# Description: Makefile de MAN/PDOC
ifeq "$(BAtail -f ../../Resultats.make.pgi13xx.ST1SE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

SHELL   = /bin/bash

CPPT    = r.gppf -chop_bang -gpp -F

DIAGNOSTIQUE = $(CURDIR)/../..
MYMAKE  = $(MAKE) DIAGNOSTIQUE=$(DIAGNOSTIQUE)

FFLAGS  =  

ENTETE  = 32
NP      = 1
OPTIL   = 1
STD     = 98

DIAGS  = $(DIAGNOSTIQUE)
LSPGM  = $(DIAGS)/src/lspgm
LSSUB  = $(DIAGS)/src/lssub
CDF2C  = $(DIAGS)/src/cdf2ccc

FTNT   = $(CPPT) -DNO_SOURCE -DAUTO_DOC

DIAG_VERSION = 6.4.2

HOSTWEB =
DIAGWEB =

default:  document

%.txt : $(LSPGM)/%.ptn
	cd $(LSPGM) ; $(MYMAKE) $@
	if ( test "$(HOSTWEB):$(DIAGWEB)" != ":" ) then \
	rsync -av $@ $(HOSTWEB):$(DIAGWEB) ; fi

%.txt : $(LSPGM)/%.ptn90
	cd  $(LSPGM) ; $(MYMAKE) $@
	if ( test "$(HOSTWEB):$(DIAGWEB)" != ":" ) then \
	rsync -av $@ $(HOSTWEB):$(DIAGWEB) ; fi

#LISTE DE TOUS LES FICHIERS .TXT

TXTDECKS = \
	addf.txt	addtext.txt	adjust.txt	anomali.txt	barx.txt	bary.txt \
	beta.txt	betao.txt	binach.txt	binning.txt	bins.txt \
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
	zxprint.txt	zxpsi.txt	zxtrig.txt	zxwave.txt	zxmav.txt \
	lssub.txt	driver.txt	cdf2ccc.txt

#LISTE DES CIBLES DU FICHIER recettes

lssub.txt	:	$(LSSUB)/lssub.ftn
	cd  $(LSSUB) ; $(MYMAKE) $@
	if ( test "$(HOSTWEB):$(DIAGWEB)" != ":" ) then \
	rsync -av $@ $(HOSTWEB):$(DIAGWEB) ; fi

cdf2ccc.txt	:	$(CDF2C)/cdf2ccc.F90
	cd  $(CDF2C) ; $(MYMAKE) $@
	if ( test "$(HOSTWEB):$(DIAGWEB)" != ":" ) then \
	rsync -av $@ $(HOSTWEB):$(DIAGWEB) ; fi

document:	info.lspgm	$(TXTDECKS)

info.lspgm	:	$(LSPGM)/driver.ptn
	for i in a b c d e f g; do \
	$(CPPT) -DNO_SOURCE -DMORE_INFO $(LSPGM)/driver.ptn | grep "^[*]$$i" | sed "s/^[*]$$i/ /" > lspgm_$$i.txt; \
	/bin/cat lspgm_$$i.txt | info2html > lspgm_$$i.html ; \
	if [[ "$(HOSTWEB):$(DIAGWEB)" != ":" ]] ; then \
	rsync -av lspgm_$$i.txt lspgm_$$i.html $(HOSTWEB):$(DIAGWEB) ; fi \
	done

web_document:
	if ( test "$(HOSTWEB):$(DIAGWEB)" != ":" ) then \
	rsync -av * $(HOSTWEB):$(DIAGWEB) ; fi
