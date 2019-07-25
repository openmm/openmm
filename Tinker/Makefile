##
###################################################################
##                                                               ##
##  Makefile for Building the TINKER Molecular Modeling Package  ##
##                                                               ##
###################################################################
##
##  Invocation Options:
##
##   1. make all              Build all the TINKER executables
##   2. make install          Move the executables to BINDIR
##   3. make create_links     Create soft links in LINKDIR
##   4. make remove_links     Remove soft links from LINKDIR
##   6. make listing          Concatenate source to tinker.txt
##   5. make clean            Delete objects and executables
##
##  Modified from the standard TINKER Makefile in March 2015 by
##  Jay Ponder by build an OpenMM-accelerated dynamic executable.
##  Requires TINKER 7.1.2 or later, and assumes an AMOEBA-capable
##  version of OpenMM is already built and present
##
##  Original version of this file is due to Peter Happersberger
##  and Jochen Buehler of the University of Konstanz, January 1998.
##  Modifications by Reece Hart & Jay Ponder, Washington University
##
##  Building TINKER from source requires libraries from the FFTW
##  Fast Fourier Transform package and (optionally) APBS Poisson-
##  Boltzmann solver. These libraries are assumed to already be
##  present, and are not built by this Makefile
##
###################################################################

###################################################################
##  Master Directory Locations; Change as Needed for Local Site  ##
###################################################################

##  TINKERDIR    TINKER Distribution Directory
##  BINDIR       Directory with TINKER Executables
##  LINKDIR      Linked Copies of TINKER Executables
##  FFTWDIR      FFTW Top-Level Directory
##  OPENMMDIR    OpenMM Distribution Directory

TINKERDIR =/home/liuchw/Softwares/tinkers
BINDIR = $(TINKERDIR)/bin
LINKDIR = /usr/local/bin

FFTWDIR =/home/liuchw/Softwares/tinkers/tinker8.2-intel13/fftw
FFTW_LIB_DIR = -L$(FFTWDIR)/lib
FFTW_LIBS = -lfftw3_threads -lfftw3

OPENMMDIR =/home/liuchw/Softwares/tinkers/Tinker-OpenMM/AMOEBAplus/bin2019
OPENMM_LIB_DIR = -L$(OPENMMDIR)/lib
OPENMM_INCLUDE_DIR = $(OPENMMDIR)/include
OPENMM_LIBS = -lOpenMM -lOpenMMAmoeba

CUDA_DIR = /usr/local/cuda-8.0
CUDA_LIB = $(CUDA_DIR)/lib64
CUDA_INCLUDE = $(CUDA_DIR)/include

NVML_INCLUDE = /usr/include/nvidia/gdk
NVML_LIB = /usr/src/gdk/nvml/lib/

CC = gcc
CXX = g++
NVCC = /usr/local/cuda/bin/nvcc
F77 = gfortran
LIBS = -lstdc++ -L. -L$(NVML_LIB) -lnvidia-ml -lfftw3_threads -lfftw3 -L$(CUDA_LIB) -lcudart -lcuda
CFLAGS = -c
F77FLAGS = -c
OPTFLAGS_C =  -O3 -ffast-math -fopenmp
OPTFLAGS_F = -O3 -ffast-math -fopenmp
CUDAFLAGS = -O3 -I $(CUDA_INCLUDE) -I$(NVML_INCLUDE)
LIBFLAGS = -crusv
RANLIB = ranlib
LINKFLAGS = $(OPTFLAGS_C) -static-libgcc 

#p################################################################
##  Should not be Necessary to Change Things Below this Point  ##
#################################################################

OBJS =	action.o \
	active.o \
	alchemy.o \
	align.o \
	analysis.o \
	analyz.o \
	analyze.o \
	angang.o \
	angbnd.o \
	angles.o \
	angpot.o \
	angtor.o \
	anneal.o \
	archive.o \
	argue.o \
	ascii.o \
	atmlst.o \
	atomid.o \
	atoms.o \
	attach.o \
	bar.o \
	basefile.o \
	bath.o \
	beeman.o \
	bicubic.o \
	bitor.o \
	bitors.o \
	bndpot.o \
	bndstr.o \
	bonds.o \
	born.o \
	bound.o \
	bounds.o \
	boxes.o \
	bussi.o \
	calendar.o \
	cell.o \
	center.o \
	charge.o \
	ctran.o \
	chgpot.o \
	chkpole.o \
	chkring.o \
	chkxyz.o \
	cholesky.o \
	chrono.o \
	chunks.o \
	clock.o \
	cluster.o \
	column.o \
	command.o \
	connect.o \
	connolly.o \
	control.o \
	correlate.o \
	couple.o \
	crystal.o \
	cspline.o \
	cutoffs.o \
	deflate.o \
	delete.o \
	deriv.o \
	diagq.o \
	diffeq.o \
	diffuse.o \
	dipole.o \
	disgeo.o \
	distgeom.o \
	dma.o \
	document.o \
	domega.o \
	dynamic.o \
	dynamic_omm.o \
	eangang.o \
	eangang1.o \
	eangang2.o \
	eangang3.o \
	eangle.o \
	eangle1.o \
	eangle2.o \
	eangle3.o \
	eangtor.o \
	eangtor1.o \
	eangtor2.o \
	eangtor3.o \
	ebond.o \
	ebond1.o \
	ebond2.o \
	ebond3.o \
	ebuck.o \
	ebuck1.o \
	ebuck2.o \
	ebuck3.o \
	echarge.o \
	echarge1.o \
	echarge2.o \
	echarge3.o \
	echgdpl.o \
	echgdpl1.o \
	echgdpl2.o \
	echgdpl3.o \
	edipole.o \
	edipole1.o \
	edipole2.o \
	edipole3.o \
	egauss.o \
	egauss1.o \
	egauss2.o \
	egauss3.o \
	egeom.o \
	egeom1.o \
	egeom2.o \
	egeom3.o \
	ehal.o \
	ehal1.o \
	ehal2.o \
	ehal3.o \
	eimprop.o \
	eimprop1.o \
	eimprop2.o \
	eimprop3.o \
	eimptor.o \
	eimptor1.o \
	eimptor2.o \
	eimptor3.o \
	elj.o \
	elj1.o \
	elj2.o \
	elj3.o \
	ect0.o \
	ect1.o \
	ect3.o \
	embed.o \
	emetal.o \
	emetal1.o \
	emetal2.o \
	emetal3.o \
	emm3hb.o \
	emm3hb1.o \
	emm3hb2.o \
	emm3hb3.o \
	empole.o \
	empole1.o \
	empole3.o \
	energi.o \
	energy.o \
	eopbend.o \
	eopbend1.o \
	eopbend2.o \
	eopbend3.o \
	eopdist.o \
	eopdist1.o \
	eopdist2.o \
	eopdist3.o \
	epitors.o \
	epitors1.o \
	epitors2.o \
	epitors3.o \
	epolar.o \
	epolar1.o \
	epolar3.o \
	erf.o \
	erxnfld.o \
	erxnfld1.o \
	erxnfld2.o \
	erxnfld3.o \
	esolv.o \
	esolv1.o \
	esolv2.o \
	esolv3.o \
	estrbnd.o \
	estrbnd1.o \
	estrbnd2.o \
	estrbnd3.o \
	estrtor.o \
	estrtor1.o \
	estrtor2.o \
	estrtor3.o \
	etors.o \
	etors1.o \
	etors2.o \
	etors3.o \
	etortor.o \
	etortor1.o \
	etortor2.o \
	etortor3.o \
	eurey.o \
	eurey1.o \
	eurey2.o \
	eurey3.o \
	evcorr.o \
	ewald.o \
	extra.o \
	extra1.o \
	extra2.o \
	extra3.o \
	faces.o \
	fatal.o \
	fft.o \
	fft3d.o \
	fftpack.o \
	field.o \
	fields.o \
	files.o \
	final.o \
	flatten.o \
	fracs.o \
	freeunit.o \
	freeze.o \
	gda.o \
	geometry.o \
	getint.o \
	getkey.o \
	getmol.o \
	getmol2.o \
	getnumb.o \
	getpdb.o \
	getprm.o \
	getref.o \
	getstring.o \
	gettext.o \
	getword.o \
	getxyz.o \
	ghmcstep.o \
	gkstuf.o \
	gpu_cards.o \
	gradient.o \
	gradrgd.o \
	gradrot.o \
	group.o \
	groups.o \
	grpline.o \
	gyrate.o \
	hescut.o \
	hessian.o \
	hessn.o \
	hessrgd.o \
	hessrot.o \
	hpmf.o \
	hybrid.o \
	ielscf.o \
	image.o \
	impose.o \
	improp.o \
	imptor.o \
	induce.o \
	inertia.o \
	inform.o \
	initatom.o \
	initial.o \
	initprm.o \
	initres.o \
	initrot.o \
	insert.o \
	intedit.o \
	inter.o \
	intxyz.o \
	invbeta.o \
	invert.o \
	iounit.o \
	jacobi.o \
	kanang.o \
	kangang.o \
	kangle.o \
	kangs.o \
	kangtor.o \
	kantor.o \
	katom.o \
	katoms.o \
	kbond.o \
	kbonds.o \
	kcharge.o \
	kchrge.o \
	kdipol.o \
	kdipole.o \
	kewald.o \
	kextra.o \
	keys.o \
	kgeom.o \
	khbond.o \
	kimprop.o \
	kimptor.o \
	kinetic.o \
	kiprop.o \
	kitors.o \
	kmetal.o \
	kmpole.o \
	kmulti.o \
	kopbend.o \
	kopbnd.o \
	kopdist.o \
	kopdst.o \
	korbit.o \
	korbs.o \
	kpitor.o \
	kpitors.o \
	kpolar.o \
	kpolr.o \
	ksolv.o \
	kstbnd.o \
	kstrbnd.o \
	kstrtor.o \
	ksttor.o \
	ktors.o \
	ktorsn.o \
	ktortor.o \
	ktrtor.o \
	kurey.o \
	kurybr.o \
	kctr.o \
	kvdw.o \
	kvdwpr.o \
	kvdws.o \
	lattice.o \
	lbfgs.o \
	light.o \
	lights.o \
	limits.o \
	linmin.o \
	makeint.o \
	makeref.o \
	makexyz.o \
	math.o \
	maxwell.o \
	mdinit.o \
	mdrest.o \
	mdsave.o \
	mdstat.o \
	mdstuf.o \
	mechanic.o \
	merck.o \
	merge.o \
	minima.o \
	minimize.o \
	minirot.o \
	minrigid.o \
	molcul.o \
	moldyn.o \
	molecule.o \
	molxyz.o \
	moment.o \
	moments.o \
	monte.o \
	mplpot.o \
	mpole.o \
	mrecip.o \
	mutant.o \
	mutate.o \
	nblist.o \
	neigh.o \
	newton.o \
	newtrot.o \
	nextarg.o \
	nexttext.o \
	nonpol.o \
	nose.o \
	nspline.o \
	nucleic.o \
	nucleo.o \
	number.o \
	numeral.o \
	numgrad.o \
	ocvm.o \
	omega.o \
	opbend.o \
	opdist.o \
	openend.o \
	openmm.o \
	openmp.o \
	optimize.o \
	optirot.o \
	optrigid.o \
	optsave.o \
	orbital.o \
	orbits.o \
	orient.o \
	orthog.o \
	output.o \
	overlap.o \
	params.o \
	path.o \
	paths.o \
	pbstuf.o \
	pdb.o \
	pdbxyz.o \
	phipsi.o \
	picalc.o \
	piorbs.o \
	pistuf.o \
	pitors.o \
	pme.o \
	pmestuff.o \
	pmpb.o \
	polar.o \
	polarize.o \
	poledit.o \
	polgrp.o \
	polpot.o \
	polymer.o \
	potent.o \
	potential.o \
	potfit.o \
	precis.o \
	precise.o \
	pressure.o \
	prmedit.o \
	prmkey.o \
	promo.o \
	protein.o \
	prtdyn.o \
	prterr.o \
	prtint.o \
	prtmol2.o \
	prtpdb.o \
	prtprm.o \
	prtseq.o \
	prtxyz.o \
	pss.o \
	pssrigid.o \
	pssrot.o \
	ptable.o \
	qmstuf.o \
	qrfact.o \
	quatfit.o \
	radial.o \
	random.o \
	rattle.o \
	readdyn.o \
	readgau.o \
	readgdma.o \
	readint.o \
	readmol.o \
	readmol2.o \
	readpdb.o \
	readprm.o \
	readseq.o \
	readxyz.o \
	refer.o \
	replica.o \
	resdue.o \
	respa.o \
	restrn.o \
	rgddyn.o \
	rgdstep.o \
	rigid.o \
	ring.o \
	rings.o \
	rmsfit.o \
	rotbnd.o \
	rotlist.o \
	rotpole.o \
	rxnfld.o \
	rxnpot.o \
	saddle.o \
	scales.o \
	scan.o \
	sdstep.o \
	search.o \
	sequen.o \
	server.o \
	shakeup.o \
	shunt.o \
	sigmoid.o \
	simplex.o \
	sizes.o \
	sktstuff.o \
	sniffer.o \
	socket.o \
	solute.o \
	sort.o \
	spacefill.o \
	spectrum.o \
	square.o \
	stodyn.o \
	strbnd.o \
	strtor.o \
	suffix.o \
	superpose.o \
	surface.o \
	surfatom.o \
	switch.o \
	sybylxyz.o \
	syntrn.o \
	tarray.o \
	temper.o \
	testgrad.o \
	testhess.o \
	testpair.o \
	testpol.o \
	testrot.o \
	timer.o \
	timerot.o \
	titles.o \
	tncg.o \
	torphase.o \
	torpot.o \
	torque.o \
	tors.o \
	torsfit.o \
	torsions.o \
	tortor.o \
	tree.o \
	trimtext.o \
	unitcell.o \
	units.o \
	uprior.o \
	urey.o \
	urypot.o \
	usage.o \
	usolve.o \
	valence.o \
	valfit.o \
	vdw.o \
	vdwpot.o \
	verlet.o \
	version.o \
	vibbig.o \
	vibrate.o \
	vibrot.o \
	vibs.o \
	virial.o \
	volume.o \
	warp.o \
	xtalfit.o \
	xtalmin.o \
	xtals.o \
	xyzatm.o \
	xyzedit.o \
	xyzint.o \
	xyzpdb.o \
	xyzsybyl.o \
	zatom.o \
	zclose.o \
	zcoord.o

EXEFILES = dynamic_omm.x

%.o: %.f
	${F77} ${F77FLAGS} ${OPTFLAGS_F} $< -o $@ 

%.o: %.c
	${CC} ${CFLAGS} ${OPTFLAGS_C} ${INCLUDEDIR} -O3 $<

%.o: %.cu
	${NVCC} ${CFLAGS} ${INCLUDEDIR} ${CUDAFLAGS} -Wno-deprecated-gpu-targets $<

%.o: %.m
	${OBJECTIVE_CC} ${OBJECTIVE_CFLAGS} ${INCLUDEDIR} ${CUDAFLAGS} $< -o $@

%.x: %.o libtinker.a
	${F77} ${LINKFLAGS} -o $@ $^ ${LIBS}; strip $@

all:${EXEFILES}

clean:
	rm -f *.o *.mod *.x

listing:
	cat *.f *.cpp > tinker.txt

ommstuff.o: ommstuff.cpp
	${CXX} ${CFLAGS} ${OPTFLAGS_C} -I${OPENMM_INCLUDE_DIR} ommstuff.cpp

dynamic_omm.x: dynamic_omm.o ommstuff.o libtinker.a
	${F77} ${LINKFLAGS} ${FFTW_LIB_DIR} ${OPENMM_LIB_DIR} -o $@ $^ ${LIBS} ${FFTW_LIBS} ${OPENMM_LIBS}; strip $@

rename:
	mv dynamic_omm.x $(BINDIR)/dynamic_omm

remove_links:
	rm -f $(LINKDIR)/dynamic_omm

create_links:
	ln -s $(BINDIR)/dynamic_omm $(LINKDIR)/dynamic_omm

libtinker.a: ${OBJS} 
	ar ${LIBFLAGS} libtinker.a \
	action.o \
	active.o \
	align.o \
	analysis.o \
	analyz.o \
	angang.o \
	angbnd.o \
	angles.o \
	angpot.o \
	angtor.o \
	argue.o \
	ascii.o \
	atmlst.o \
	atomid.o \
	atoms.o \
	attach.o \
	basefile.o \
	bath.o \
	beeman.o \
	bicubic.o \
	bitor.o \
	bitors.o \
	bndpot.o \
	bndstr.o \
	bonds.o \
	born.o \
	bound.o \
	bounds.o \
	boxes.o \
	bussi.o \
	calendar.o \
	cell.o \
	center.o \
	charge.o \
	chgpot.o \
	chkpole.o \
	chkring.o \
	chkxyz.o \
	cholesky.o \
	chrono.o \
	chunks.o \
	clock.o \
	cluster.o \
	column.o \
	command.o \
	connect.o \
	connolly.o \
	control.o \
	couple.o \
	ctran.o \
	cspline.o \
	cutoffs.o \
	deflate.o \
	delete.o \
	deriv.o \
	diagq.o \
	diffeq.o \
	dipole.o \
	disgeo.o \
	dma.o \
	domega.o \
	eangang.o \
	eangang1.o \
	eangang2.o \
	eangang3.o \
	eangle.o \
	eangle1.o \
	eangle2.o \
	eangle3.o \
	eangtor.o \
	eangtor1.o \
	eangtor2.o \
	eangtor3.o \
	ebond.o \
	ebond1.o \
	ebond2.o \
	ebond3.o \
	ebuck.o \
	ebuck1.o \
	ebuck2.o \
	ebuck3.o \
	echarge.o \
	echarge1.o \
	echarge2.o \
	echarge3.o \
	echgdpl.o \
	echgdpl1.o \
	echgdpl2.o \
	echgdpl3.o \
	edipole.o \
	edipole1.o \
	edipole2.o \
	edipole3.o \
	egauss.o \
	egauss1.o \
	egauss2.o \
	egauss3.o \
	egeom.o \
	egeom1.o \
	egeom2.o \
	egeom3.o \
	ehal.o \
	ehal1.o \
	ehal2.o \
	ehal3.o \
	eimprop.o \
	eimprop1.o \
	eimprop2.o \
	eimprop3.o \
	eimptor.o \
	eimptor1.o \
	eimptor2.o \
	eimptor3.o \
	elj.o \
	elj1.o \
	elj2.o \
	elj3.o \
	ect0.o \
	ect1.o \
	ect3.o \
	embed.o \
	emetal.o \
	emetal1.o \
	emetal2.o \
	emetal3.o \
	emm3hb.o \
	emm3hb1.o \
	emm3hb2.o \
	emm3hb3.o \
	empole.o \
	empole1.o \
	empole3.o \
	energi.o \
	energy.o \
	eopbend.o \
	eopbend1.o \
	eopbend2.o \
	eopbend3.o \
	eopdist.o \
	eopdist1.o \
	eopdist2.o \
	eopdist3.o \
	epitors.o \
	epitors1.o \
	epitors2.o \
	epitors3.o \
	epolar.o \
	epolar1.o \
	epolar3.o \
	erf.o \
	erxnfld.o \
	erxnfld1.o \
	erxnfld2.o \
	erxnfld3.o \
	esolv.o \
	esolv1.o \
	esolv2.o \
	esolv3.o \
	estrbnd.o \
	estrbnd1.o \
	estrbnd2.o \
	estrbnd3.o \
	estrtor.o \
	estrtor1.o \
	estrtor2.o \
	estrtor3.o \
	etors.o \
	etors1.o \
	etors2.o \
	etors3.o \
	etortor.o \
	etortor1.o \
	etortor2.o \
	etortor3.o \
	eurey.o \
	eurey1.o \
	eurey2.o \
	eurey3.o \
	evcorr.o \
	ewald.o \
	extra.o \
	extra1.o \
	extra2.o \
	extra3.o \
	faces.o \
	fatal.o \
	fft.o \
	fft3d.o \
	fftpack.o \
	field.o \
	fields.o \
	files.o \
	final.o \
	flatten.o \
	fracs.o \
	freeunit.o \
	freeze.o \
	geometry.o \
	getint.o \
	getkey.o \
	getmol.o \
	getmol2.o \
	getnumb.o \
	getpdb.o \
	getprm.o \
	getref.o \
	getstring.o \
	gettext.o \
	getword.o \
	getxyz.o \
	ghmcstep.o \
	gkstuf.o \
	gpu_cards.o \
	$(GPU_USAGE_OBJS) \
	gradient.o \
	gradrgd.o \
	gradrot.o \
	group.o \
	groups.o \
	grpline.o \
	gyrate.o \
	hescut.o \
	hessian.o \
	hessn.o \
	hessrgd.o \
	hessrot.o \
	hpmf.o \
	hybrid.o \
	ielscf.o \
	image.o \
	impose.o \
	improp.o \
	imptor.o \
	induce.o \
	inertia.o \
	inform.o \
	initatom.o \
	initial.o \
	initprm.o \
	initres.o \
	initrot.o \
	insert.o \
	inter.o \
	invbeta.o \
	invert.o \
	iounit.o \
	jacobi.o \
	kanang.o \
	kangang.o \
	kangle.o \
	kangs.o \
	kangtor.o \
	kantor.o \
	katom.o \
	katoms.o \
	kbond.o \
	kbonds.o \
	kcharge.o \
	kchrge.o \
	kdipol.o \
	kdipole.o \
	kewald.o \
	kextra.o \
	keys.o \
	kgeom.o \
	khbond.o \
	kimprop.o \
	kimptor.o \
	kinetic.o \
	kiprop.o \
	kitors.o \
	kmetal.o \
	kmpole.o \
	kmulti.o \
	kopbend.o \
	kopbnd.o \
	kopdist.o \
	kopdst.o \
	korbit.o \
	korbs.o \
	kpitor.o \
	kpitors.o \
	kpolar.o \
	kpolr.o \
	ksolv.o \
	kstbnd.o \
	kstrbnd.o \
	kstrtor.o \
	ksttor.o \
	ktors.o \
	ktorsn.o \
	ktortor.o \
	ktrtor.o \
	kurey.o \
	kurybr.o \
	kctr.o \
	kvdw.o \
	kvdwpr.o \
	kvdws.o \
	lattice.o \
	lbfgs.o \
	light.o \
	lights.o \
	limits.o \
	linmin.o \
	makeint.o \
	makeref.o \
	makexyz.o \
	math.o \
	maxwell.o \
	mdinit.o \
	mdrest.o \
	mdsave.o \
	mdstat.o \
	mdstuf.o \
	mechanic.o \
	merck.o \
	merge.o \
	minima.o \
	molcul.o \
	moldyn.o \
	molecule.o \
	moment.o \
	moments.o \
	mplpot.o \
	mpole.o \
	mrecip.o \
	mutant.o \
	mutate.o \
	nblist.o \
	neigh.o \
	nextarg.o \
	nexttext.o \
	nonpol.o \
	nose.o \
	nspline.o \
	nucleo.o \
	number.o \
	numeral.o \
	numgrad.o \
	ocvm.o \
	omega.o \
	opbend.o \
	opdist.o \
	openend.o \
	openmm.o \
	openmp.o \
	optsave.o \
	orbital.o \
	orbits.o \
	orient.o \
	orthog.o \
	output.o \
	overlap.o \
	params.o \
	paths.o \
	pbstuf.o \
	pdb.o \
	phipsi.o \
	picalc.o \
	piorbs.o \
	pistuf.o \
	pitors.o \
	pme.o \
	pmestuff.o \
	pmpb.o \
	polar.o \
	polgrp.o \
	polpot.o \
	polymer.o \
	potent.o \
	potfit.o \
	precis.o \
	precise.o \
	pressure.o \
	prmkey.o \
	promo.o \
	prtdyn.o \
	prterr.o \
	prtint.o \
	prtmol2.o \
	prtpdb.o \
	prtprm.o \
	prtseq.o \
	prtxyz.o \
	ptable.o \
	qmstuf.o \
	qrfact.o \
	quatfit.o \
	random.o \
	rattle.o \
	readdyn.o \
	readgau.o \
	readgdma.o \
	readint.o \
	readmol.o \
	readmol2.o \
	readpdb.o \
	readprm.o \
	readseq.o \
	readxyz.o \
	refer.o \
	replica.o \
	resdue.o \
	respa.o \
	restrn.o \
	rgddyn.o \
	rgdstep.o \
	rigid.o \
	ring.o \
	rings.o \
	rmsfit.o \
	rotbnd.o \
	rotlist.o \
	rotpole.o \
	rxnfld.o \
	rxnpot.o \
	scales.o \
	sdstep.o \
	search.o \
	sequen.o \
	server.o \
	shakeup.o \
	shunt.o \
	sigmoid.o \
	simplex.o \
	sizes.o \
	sktstuff.o \
	socket.o \
	solute.o \
	sort.o \
	square.o \
	stodyn.o \
	strbnd.o \
	strtor.o \
	suffix.o \
	surface.o \
	surfatom.o \
	switch.o \
	syntrn.o \
	tarray.o \
	temper.o \
	titles.o \
	tncg.o \
	torphase.o \
	torpot.o \
	torque.o \
	tors.o \
	torsions.o \
	tortor.o \
	tree.o \
	trimtext.o \
	unitcell.o \
	units.o \
	uprior.o \
	urey.o \
	urypot.o \
	usage.o \
	usolve.o \
	valfit.o \
	vdw.o \
	vdwpot.o \
	verlet.o \
	version.o \
	vibs.o \
	virial.o \
	volume.o \
	warp.o \
	xtals.o \
	xyzatm.o \
	zatom.o \
	zclose.o \
	zcoord.o
	${RANLIB} libtinker.a

###############################################################
##  Next Section has Explicit Dependencies on Include Files  ##
###############################################################

dynamic_omm.o: angbnd.o angpot.o angtor.o atomid.o atoms.o bath.o bitor.o bndpot.o bndstr.o bound.o boxes.o cell.o charge.o chgpot.o couple.o ctran.o deriv.o energi.o ewald.o freeze.o group.o imptor.o inform.o iounit.o keys.o ktrtor.o kvdws.o limits.o mdstuf.o molcul.o moldyn.o mplpot.o mpole.o mutant.o nonpol.o opbend.o openmm.o openmp.o pitors.o pme.o polar.o polgrp.o polpot.o potent.o restrn.o sizes.o solute.o stodyn.o strbnd.o strtor.o torpot.o tors.o tortor.o units.o urey.o urypot.o usage.o vdw.o vdwpot.o
gpu_cards.o: $(GPU_USAGE_OBJS)
openmm.o:

action.o:
active.o: atoms.o inform.o iounit.o keys.o sizes.o usage.o
alchemy.o: analyz.o atoms.o energi.o files.o inform.o iounit.o katoms.o mutant.o potent.o sizes.o units.o usage.o
align.o:
analysis.o: analyz.o atoms.o bound.o energi.o group.o inter.o iounit.o limits.o potent.o sizes.o vdwpot.o
analyz.o:
analyze.o: action.o analyz.o angang.o angbnd.o angpot.o angtor.o atomid.o atoms.o bitor.o bndstr.o bound.o boxes.o charge.o chgpot.o couple.o ctran.o dipole.o energi.o ewald.o fields.o files.o improp.o imptor.o inform.o inter.o iounit.o korbs.o ktrtor.o kvdws.o limits.o math.o molcul.o moment.o mpole.o opbend.o opdist.o piorbs.o pistuf.o pitors.o pme.o polar.o polgrp.o potent.o sizes.o solute.o strbnd.o strtor.o tors.o tortor.o units.o urey.o vdw.o vdwpot.o virial.o
angang.o:
angbnd.o:
angles.o: angbnd.o atmlst.o atoms.o couple.o iounit.o sizes.o
angpot.o:
angtor.o:
anneal.o: atomid.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o mdstuf.o potent.o sizes.o solute.o usage.o warp.o
archive.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o iounit.o sizes.o titles.o usage.o
argue.o:
ascii.o:
atmlst.o:
atomid.o: sizes.o
atoms.o: sizes.o
attach.o: atoms.o couple.o iounit.o sizes.o
bar.o: boxes.o files.o inform.o iounit.o keys.o titles.o units.o
basefile.o: ascii.o files.o
bath.o:
beeman.o: atomid.o atoms.o freeze.o ielscf.o mdstuf.o moldyn.o polar.o sizes.o units.o usage.o
bicubic.o:
bitor.o:
bitors.o: angbnd.o atoms.o bitor.o couple.o iounit.o sizes.o
bndpot.o:
bndstr.o:
bonds.o: atmlst.o atoms.o bndstr.o couple.o iounit.o sizes.o
born.o: atomid.o atoms.o bath.o chgpot.o couple.o deriv.o inform.o iounit.o math.o pbstuf.o sizes.o solute.o virial.o
bound.o:
bounds.o: atomid.o atoms.o boxes.o molcul.o sizes.o
boxes.o:
bussi.o: atomid.o atoms.o bath.o boxes.o freeze.o ielscf.o mdstuf.o moldyn.o polar.o sizes.o units.o usage.o
calendar.o:
cell.o:
center.o: align.o
charge.o:
ctran.o: sizes.o
chgpot.o:
chkpole.o: atoms.o mpole.o sizes.o
chkring.o: couple.o
chkxyz.o: atoms.o iounit.o sizes.o
cholesky.o:
chrono.o:
chunks.o:
clock.o: chrono.o
cluster.o: atomid.o atoms.o bound.o group.o inform.o iounit.o keys.o limits.o molcul.o sizes.o
column.o:
command.o: argue.o
connect.o: atoms.o couple.o sizes.o zclose.o zcoord.o
connolly.o: atoms.o faces.o inform.o iounit.o math.o sizes.o
control.o: argue.o inform.o keys.o output.o
correlate.o: ascii.o atomid.o atoms.o files.o inform.o iounit.o sizes.o
couple.o: sizes.o
crystal.o: atomid.o atoms.o bound.o boxes.o couple.o files.o iounit.o math.o molcul.o sizes.o
cspline.o: iounit.o
cutoffs.o: atoms.o bound.o hescut.o keys.o limits.o neigh.o polpot.o sizes.o tarray.o
deflate.o: iounit.o
delete.o: atomid.o atoms.o couple.o inform.o iounit.o sizes.o
deriv.o:
diagq.o:
diffeq.o: atoms.o iounit.o math.o sizes.o warp.o
diffuse.o: atomid.o atoms.o bound.o inform.o iounit.o molcul.o sizes.o usage.o
dipole.o:
disgeo.o:
distgeom.o: angbnd.o atomid.o atoms.o bndstr.o couple.o disgeo.o files.o inform.o iounit.o kvdws.o math.o refer.o restrn.o sizes.o tors.o
dma.o:
document.o: iounit.o
domega.o:
dynamic.o: atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o sizes.o solute.o stodyn.o usage.o
eangang.o: angang.o angbnd.o angpot.o atoms.o bound.o energi.o group.o math.o sizes.o usage.o
eangang1.o: angang.o angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o sizes.o usage.o virial.o
eangang2.o: angang.o angbnd.o angpot.o atoms.o bound.o group.o hessn.o math.o sizes.o
eangang3.o: action.o analyz.o angang.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o usage.o
eangle.o: angbnd.o angpot.o atoms.o bound.o energi.o group.o math.o sizes.o usage.o
eangle1.o: angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o sizes.o usage.o virial.o
eangle2.o: angbnd.o angpot.o atoms.o bound.o group.o hessn.o math.o sizes.o
eangle3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o usage.o
eangtor.o: angbnd.o angtor.o atoms.o bound.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o
eangtor1.o: angbnd.o angtor.o atoms.o bound.o deriv.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o virial.o
eangtor2.o: angbnd.o angtor.o atoms.o bound.o group.o hessn.o math.o sizes.o torpot.o tors.o
eangtor3.o: action.o analyz.o angbnd.o angtor.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o torpot.o tors.o usage.o
ebond.o: atoms.o bndpot.o bndstr.o bound.o energi.o group.o sizes.o usage.o
ebond1.o: atoms.o bndpot.o bndstr.o bound.o deriv.o energi.o group.o sizes.o usage.o virial.o
ebond2.o: atmlst.o atoms.o bndpot.o bndstr.o bound.o couple.o group.o hessn.o sizes.o
ebond3.o: action.o analyz.o atomid.o atoms.o bndpot.o bndstr.o bound.o energi.o group.o inform.o iounit.o sizes.o usage.o
ebuck.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o iounit.o light.o limits.o math.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
ebuck1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o warp.o
ebuck2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o iounit.o math.o shunt.o sizes.o vdw.o vdwpot.o warp.o
ebuck3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
echarge.o: atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o energi.o ewald.o group.o iounit.o light.o limits.o math.o neigh.o pme.o shunt.o sizes.o usage.o warp.o
echarge1.o: atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o light.o limits.o math.o molcul.o neigh.o pme.o shunt.o sizes.o usage.o virial.o warp.o
echarge2.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o ewald.o group.o hessn.o limits.o math.o shunt.o sizes.o warp.o
echarge3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o energi.o ewald.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o warp.o
echgdpl.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o energi.o group.o shunt.o sizes.o units.o usage.o
echgdpl1.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o deriv.o dipole.o energi.o group.o inter.o molcul.o shunt.o sizes.o units.o usage.o virial.o
echgdpl2.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o group.o hessn.o shunt.o sizes.o units.o
echgdpl3.o: action.o analyz.o atomid.o atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o energi.o group.o inform.o inter.o iounit.o molcul.o shunt.o sizes.o units.o usage.o
edipole.o: atoms.o bound.o cell.o chgpot.o dipole.o energi.o group.o shunt.o sizes.o units.o usage.o
edipole1.o: atoms.o bound.o cell.o chgpot.o deriv.o dipole.o energi.o group.o inter.o molcul.o shunt.o sizes.o units.o usage.o virial.o
edipole2.o: atoms.o bound.o cell.o chgpot.o dipole.o group.o hessn.o shunt.o sizes.o units.o
edipole3.o: action.o analyz.o atomid.o atoms.o bound.o cell.o chgpot.o dipole.o energi.o group.o inform.o inter.o iounit.o molcul.o shunt.o sizes.o units.o usage.o
egauss.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o math.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
egauss1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o warp.o
egauss2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o shunt.o sizes.o vdw.o vdwpot.o warp.o
egauss3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
egeom.o: atomid.o atoms.o bound.o energi.o group.o math.o molcul.o restrn.o sizes.o usage.o
egeom1.o: atomid.o atoms.o bound.o deriv.o energi.o group.o inter.o math.o molcul.o restrn.o sizes.o usage.o virial.o
egeom2.o: atomid.o atoms.o bound.o deriv.o group.o hessn.o math.o molcul.o restrn.o sizes.o
egeom3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o inter.o iounit.o math.o molcul.o restrn.o sizes.o usage.o
ehal.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
ehal1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o molcul.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o
ehal2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o shunt.o sizes.o vdw.o vdwpot.o
ehal3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o molcul.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
eimprop.o: atoms.o bound.o energi.o group.o improp.o math.o sizes.o torpot.o usage.o
eimprop1.o: atoms.o bound.o deriv.o energi.o group.o improp.o math.o sizes.o torpot.o usage.o virial.o
eimprop2.o: atoms.o bound.o group.o hessn.o improp.o math.o sizes.o torpot.o
eimprop3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o improp.o inform.o iounit.o math.o sizes.o torpot.o usage.o
eimptor.o: atoms.o bound.o energi.o group.o imptor.o sizes.o torpot.o usage.o
eimptor1.o: atoms.o bound.o deriv.o energi.o group.o imptor.o sizes.o torpot.o usage.o virial.o
eimptor2.o: atoms.o bound.o group.o hessn.o imptor.o sizes.o torpot.o
eimptor3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o imptor.o inform.o iounit.o math.o sizes.o torpot.o usage.o
elj.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o math.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
elj1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o warp.o
elj2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o math.o shunt.o sizes.o vdw.o vdwpot.o warp.o
elj3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
ect0.o: atomid.o atoms.o bound.o sizes.o ctran.o cell.o couple.o energi.o group.o shunt.o usage.o 
ect1.o: sizes.o limits.o warp.o ctran.o atomid.o atoms.o bound.o cell.o couple.o energi.o group.o shunt.o usage.o virial.o deriv.o inter.o molcul.o  
ect3.o: action.o analyz.o sizes.o limits.o warp.o ctran.o atomid.o atoms.o bound.o cell.o couple.o energi.o group.o molcul.o inter.o iounit.o inform.o shunt.o usage.o 
embed.o: angbnd.o atoms.o bndstr.o couple.o disgeo.o files.o inform.o iounit.o keys.o light.o math.o minima.o output.o refer.o restrn.o sizes.o tors.o units.o
emetal.o: atomid.o atoms.o couple.o energi.o kchrge.o sizes.o
emetal1.o: atomid.o atoms.o couple.o deriv.o energi.o kchrge.o sizes.o
emetal2.o:
emetal3.o: action.o analyz.o atomid.o atoms.o energi.o kchrge.o sizes.o
emm3hb.o: atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o energi.o group.o light.o limits.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
emm3hb1.o: atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o
emm3hb2.o: atmlst.o atomid.o atoms.o bndstr.o bound.o cell.o chgpot.o couple.o group.o hessn.o shunt.o sizes.o vdw.o vdwpot.o
emm3hb3.o: action.o analyz.o atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
empole.o: atoms.o bound.o boxes.o cell.o chgpot.o ctran.o couple.o energi.o ewald.o group.o limits.o math.o mplpot.o mpole.o mrecip.o neigh.o pme.o potent.o shunt.o sizes.o usage.o
empole1.o: atoms.o bound.o boxes.o cell.o chgpot.o ctran.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o mrecip.o neigh.o pme.o shunt.o sizes.o usage.o virial.o
empole3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o chgpot.o ctran.o couple.o energi.o ewald.o group.o inform.o inter.o iounit.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o potent.o shunt.o sizes.o usage.o
energi.o:
energy.o: bound.o energi.o iounit.o limits.o potent.o rigid.o sizes.o vdwpot.o
eopbend.o: angbnd.o angpot.o atoms.o bound.o energi.o fields.o group.o math.o opbend.o sizes.o usage.o
eopbend1.o: angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o opbend.o sizes.o usage.o virial.o
eopbend2.o: angbnd.o angpot.o atoms.o bound.o group.o hessn.o math.o opbend.o sizes.o
eopbend3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o opbend.o sizes.o usage.o
eopdist.o: angpot.o atoms.o bound.o energi.o group.o opdist.o sizes.o usage.o
eopdist1.o: angpot.o atoms.o bound.o deriv.o energi.o group.o opdist.o sizes.o usage.o virial.o
eopdist2.o: angpot.o atoms.o bound.o group.o hessn.o opdist.o sizes.o usage.o
eopdist3.o: action.o analyz.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o opdist.o sizes.o usage.o
epitors.o: atoms.o bound.o energi.o group.o pitors.o sizes.o torpot.o usage.o
epitors1.o: atoms.o bound.o deriv.o energi.o group.o pitors.o sizes.o torpot.o usage.o virial.o
epitors2.o: angbnd.o atoms.o bound.o deriv.o group.o hessn.o pitors.o sizes.o torpot.o usage.o
epitors3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o pitors.o sizes.o torpot.o usage.o
epolar.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o energi.o ewald.o limits.o math.o mpole.o mrecip.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o units.o
epolar1.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o inter.o limits.o math.o molcul.o mpole.o mrecip.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o virial.o
epolar3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o chgpot.o couple.o energi.o ewald.o inform.o inter.o iounit.o limits.o math.o molcul.o mpole.o neigh.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o units.o
erf.o: iounit.o math.o
erxnfld.o: atoms.o chgpot.o energi.o mpole.o rxnfld.o rxnpot.o shunt.o sizes.o usage.o
erxnfld1.o: atoms.o deriv.o energi.o sizes.o
erxnfld2.o:
erxnfld3.o: action.o analyz.o atomid.o atoms.o chgpot.o energi.o inform.o iounit.o mpole.o shunt.o sizes.o usage.o
esolv.o: atomid.o atoms.o bound.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o kvdws.o limits.o math.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o usage.o vdw.o warp.o
esolv1.o: atomid.o atoms.o bound.o boxes.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o inter.o iounit.o kvdws.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o usage.o vdw.o virial.o warp.o
esolv2.o: atoms.o charge.o chgpot.o deriv.o hessn.o math.o mpole.o potent.o shunt.o sizes.o solute.o warp.o
esolv3.o: action.o analyz.o atomid.o atoms.o bound.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o inform.o inter.o iounit.o kvdws.o limits.o math.o molcul.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o usage.o vdw.o warp.o
estrbnd.o: angbnd.o angpot.o atoms.o bndstr.o bound.o energi.o group.o math.o sizes.o strbnd.o usage.o
estrbnd1.o: angbnd.o angpot.o atoms.o bndstr.o bound.o deriv.o energi.o group.o math.o sizes.o strbnd.o usage.o virial.o
estrbnd2.o: angbnd.o angpot.o atoms.o bndstr.o bound.o group.o hessn.o math.o sizes.o strbnd.o
estrbnd3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bndstr.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o strbnd.o usage.o
estrtor.o: atoms.o bndstr.o bound.o energi.o group.o sizes.o strtor.o torpot.o tors.o usage.o
estrtor1.o: atoms.o bndstr.o bound.o deriv.o energi.o group.o sizes.o strtor.o torpot.o tors.o usage.o virial.o
estrtor2.o: atoms.o bndstr.o bound.o group.o hessn.o sizes.o strtor.o torpot.o tors.o
estrtor3.o: action.o analyz.o atomid.o atoms.o bndstr.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o strtor.o torpot.o tors.o usage.o
etors.o: atoms.o bound.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o warp.o
etors1.o: atoms.o bound.o deriv.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o virial.o warp.o
etors2.o: atoms.o bound.o group.o hessn.o math.o sizes.o torpot.o tors.o warp.o
etors3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o torpot.o tors.o usage.o warp.o
etortor.o: atomid.o atoms.o bitor.o bound.o couple.o energi.o group.o ktrtor.o math.o sizes.o torpot.o tortor.o usage.o
etortor1.o: atoms.o bitor.o bound.o deriv.o energi.o group.o ktrtor.o math.o sizes.o torpot.o tortor.o usage.o virial.o
etortor2.o: atoms.o bitor.o bound.o group.o hessn.o ktrtor.o math.o sizes.o torpot.o tortor.o units.o
etortor3.o: action.o analyz.o atoms.o bitor.o bound.o energi.o group.o inform.o iounit.o ktrtor.o math.o sizes.o torpot.o tortor.o usage.o
eurey.o: atoms.o bound.o energi.o group.o sizes.o urey.o urypot.o usage.o
eurey1.o: atoms.o bound.o deriv.o energi.o group.o sizes.o urey.o urypot.o usage.o virial.o
eurey2.o: atoms.o bound.o couple.o group.o hessn.o sizes.o urey.o urypot.o
eurey3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o sizes.o urey.o urypot.o usage.o
evcorr.o: atoms.o bound.o boxes.o limits.o math.o mutant.o shunt.o sizes.o vdw.o vdwpot.o
ewald.o:
extra.o: energi.o
extra1.o: atoms.o deriv.o energi.o sizes.o
extra2.o: atoms.o hessn.o sizes.o
extra3.o: action.o analyz.o atoms.o energi.o sizes.o
faces.o:
fatal.o: iounit.o
fft.o:
fft3d.o: fft.o openmp.o pme.o
fftpack.o: math.o
field.o: keys.o potent.o
fields.o:
files.o:
final.o: align.o analyz.o angang.o angbnd.o angtor.o atmlst.o bitor.o bndstr.o cell.o charge.o chunks.o couple.o ctran.o deriv.o dipole.o disgeo.o domega.o faces.o fracs.o freeze.o group.o hessn.o hpmf.o ielscf.o improp.o imptor.o inform.o iounit.o katoms.o kchrge.o kpolr.o kvdws.o light.o merck.o molcul.o moldyn.o mpole.o mrecip.o mutant.o neigh.o nonpol.o omega.o opbend.o opdist.o orbits.o paths.o pbstuf.o pdb.o piorbs.o pistuf.o pitors.o pme.o polar.o polgrp.o potfit.o qmstuf.o refer.o restrn.o rgddyn.o rigid.o ring.o rotbnd.o sizes.o socket.o solute.o stodyn.o strbnd.o strtor.o syntrn.o tarray.o tors.o tortor.o uprior.o urey.o usage.o usolve.o vdw.o vibs.o warp.o
flatten.o: atoms.o fields.o inform.o iounit.o keys.o sizes.o warp.o
fracs.o:
freeunit.o: iounit.o
freeze.o:
gda.o: atoms.o files.o iounit.o minima.o potent.o sizes.o vdwpot.o warp.o
geometry.o: atoms.o math.o sizes.o
getint.o: atoms.o inform.o iounit.o output.o sizes.o
getkey.o: argue.o files.o iounit.o keys.o openmp.o
getmol.o: files.o iounit.o
getmol2.o: files.o iounit.o
getnumb.o: ascii.o
getpdb.o: iounit.o
getprm.o: files.o iounit.o keys.o params.o
getref.o: atomid.o atoms.o boxes.o couple.o files.o refer.o sizes.o titles.o
getstring.o: ascii.o
gettext.o: ascii.o
getword.o: ascii.o
getxyz.o: inform.o iounit.o output.o
ghmcstep.o: atomid.o atoms.o bath.o freeze.o iounit.o mdstuf.o moldyn.o sizes.o stodyn.o units.o usage.o virial.o
gkstuf.o: sizes.o
gradient.o: atoms.o bound.o couple.o deriv.o energi.o inter.o iounit.o limits.o potent.o rigid.o sizes.o vdwpot.o virial.o
gradrgd.o: atoms.o group.o rigid.o sizes.o
gradrot.o: atoms.o deriv.o domega.o omega.o potent.o rotbnd.o sizes.o
group.o:
groups.o: group.o
grpline.o: atomid.o atoms.o group.o rgddyn.o sizes.o
gyrate.o: atoms.o sizes.o usage.o
hescut.o:
hessian.o: atoms.o bound.o couple.o hescut.o hessn.o inform.o iounit.o limits.o mpole.o potent.o rigid.o sizes.o usage.o vdw.o vdwpot.o
hessn.o:
hessrgd.o: atoms.o group.o rigid.o sizes.o
hessrot.o: math.o omega.o sizes.o zcoord.o
hpmf.o:
hybrid.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o charge.o couple.o dipole.o imptor.o inform.o iounit.o kangs.o katoms.o kbonds.o kchrge.o kdipol.o kitors.o kstbnd.o ksttor.o ktorsn.o kvdws.o math.o mutant.o sizes.o strbnd.o strtor.o tors.o vdw.o vdwpot.o
ielscf.o:
image.o: boxes.o cell.o sizes.o
impose.o: align.o inform.o iounit.o sizes.o
improp.o:
imptor.o:
induce.o: atoms.o bound.o boxes.o cell.o couple.o ewald.o gkstuf.o group.o ielscf.o inform.o iounit.o limits.o math.o mpole.o neigh.o openmp.o pbstuf.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o tarray.o units.o uprior.o usolve.o
inertia.o: atomid.o atoms.o iounit.o math.o sizes.o
inform.o:
initatom.o: ptable.o
initial.o: align.o atoms.o bath.o bound.o boxes.o cell.o files.o group.o inform.o iounit.o keys.o linmin.o minima.o molcul.o mutant.o neigh.o openmp.o output.o params.o pdb.o precis.o rigid.o scales.o sequen.o sizes.o socket.o virial.o warp.o zclose.o
initprm.o: angpot.o bndpot.o chgpot.o ctran.o fields.o kanang.o kangs.o kantor.o katoms.o kbonds.o kchrge.o kdipol.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o math.o merck.o mplpot.o polpot.o rxnpot.o sizes.o solute.o torpot.o units.o urypot.o vdwpot.o
initres.o: resdue.o
initrot.o: atoms.o couple.o group.o inform.o iounit.o math.o omega.o potent.o restrn.o rotbnd.o sizes.o usage.o zcoord.o
insert.o: atomid.o atoms.o couple.o inform.o iounit.o sizes.o
intedit.o: atomid.o atoms.o files.o iounit.o katoms.o sizes.o zcoord.o
inter.o:
intxyz.o: files.o iounit.o titles.o
invbeta.o:
invert.o: iounit.o
iounit.o:
jacobi.o: iounit.o
kanang.o: sizes.o
kangang.o: angang.o angbnd.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o kanang.o keys.o potent.o sizes.o tors.o
kangle.o: angbnd.o angpot.o atomid.o atoms.o bndstr.o couple.o fields.o inform.o iounit.o kangs.o keys.o merck.o potent.o ring.o sizes.o usage.o
kangs.o:
kangtor.o: angtor.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o kantor.o keys.o potent.o sizes.o tors.o
kantor.o:
katom.o: atomid.o atoms.o couple.o inform.o iounit.o katoms.o keys.o sizes.o
katoms.o:
kbond.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o couple.o fields.o inform.o iounit.o kbonds.o keys.o merck.o potent.o sizes.o tors.o usage.o
kbonds.o:
kcharge.o: atomid.o atoms.o charge.o chgpot.o couple.o fields.o inform.o iounit.o kchrge.o keys.o merck.o potent.o sizes.o
kchrge.o:
kdipol.o:
kdipole.o: atmlst.o atoms.o bndstr.o couple.o dipole.o inform.o iounit.o kdipol.o keys.o potent.o sizes.o
kewald.o: atoms.o bound.o boxes.o chunks.o ewald.o fft.o inform.o iounit.o keys.o limits.o math.o openmp.o pme.o sizes.o
kextra.o: 
keys.o:
kgeom.o: atomid.o atoms.o bound.o couple.o group.o iounit.o keys.o molcul.o potent.o restrn.o sizes.o
khbond.o:
kimprop.o: atomid.o atoms.o couple.o improp.o inform.o iounit.o keys.o kiprop.o potent.o sizes.o tors.o
kimptor.o: atomid.o atoms.o couple.o imptor.o inform.o iounit.o keys.o kitors.o math.o potent.o sizes.o tors.o
kinetic.o: atomid.o atoms.o bath.o group.o ielscf.o mdstuf.o moldyn.o rgddyn.o sizes.o units.o usage.o
kiprop.o:
kitors.o:
kmetal.o:
kmpole.o: atoms.o couple.o inform.o iounit.o keys.o kmulti.o mpole.o polar.o polgrp.o potent.o sizes.o units.o
kmulti.o:
kopbend.o: angbnd.o angpot.o atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o kopbnd.o merck.o opbend.o potent.o sizes.o usage.o
kopbnd.o:
kopdist.o: angbnd.o angpot.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o keys.o kopdst.o opdist.o potent.o sizes.o
kopdst.o:
korbit.o: atomid.o atoms.o bndstr.o inform.o iounit.o keys.o korbs.o orbits.o piorbs.o pistuf.o sizes.o tors.o units.o
korbs.o: sizes.o
kpitor.o:
kpitors.o: atomid.o atoms.o bndstr.o couple.o inform.o iounit.o keys.o kpitor.o pitors.o potent.o sizes.o tors.o
kpolar.o: atoms.o couple.o inform.o iounit.o keys.o kpolr.o mpole.o neigh.o polar.o polgrp.o polpot.o potent.o sizes.o usolve.o
kpolr.o:
ksolv.o: angbnd.o atmlst.o atomid.o atoms.o bath.o bndstr.o chgpot.o couple.o gkstuf.o hpmf.o inform.o iounit.o keys.o kvdws.o math.o nonpol.o pbstuf.o polar.o potent.o ptable.o sizes.o solute.o
kstbnd.o:
kstrbnd.o: angbnd.o angpot.o atmlst.o atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o kstbnd.o merck.o potent.o ring.o sizes.o strbnd.o
kstrtor.o: atmlst.o atomid.o atoms.o couple.o inform.o iounit.o keys.o ksttor.o potent.o sizes.o strtor.o tors.o
ksttor.o:
ktors.o: atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o ktorsn.o math.o merck.o potent.o ring.o sizes.o tors.o usage.o
ktorsn.o:
ktortor.o: atomid.o atoms.o bitor.o inform.o iounit.o keys.o ktrtor.o potent.o sizes.o tortor.o
ktrtor.o:
kurey.o: angbnd.o atomid.o atoms.o inform.o iounit.o keys.o kurybr.o potent.o sizes.o urey.o
kurybr.o:
kctr.o: atomid.o atoms.o couple.o ctran.o fields.o inform.o iounit.o keys.o math.o potent.o sizes.o 
kvdw.o: atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o khbond.o kvdwpr.o kvdws.o math.o merck.o potent.o sizes.o vdw.o vdwpot.o
kvdwpr.o:
kvdws.o:
lattice.o: boxes.o cell.o inform.o iounit.o math.o
lbfgs.o: inform.o iounit.o keys.o linmin.o math.o minima.o output.o scales.o
light.o:
lights.o: bound.o boxes.o cell.o iounit.o light.o sizes.o
limits.o:
linmin.o:
makeint.o: atoms.o couple.o inform.o iounit.o math.o sizes.o zclose.o zcoord.o
makeref.o: atomid.o atoms.o boxes.o couple.o files.o refer.o sizes.o titles.o
makexyz.o: atoms.o sizes.o zcoord.o
math.o:
maxwell.o: units.o
mdinit.o: atomid.o atoms.o bath.o bound.o couple.o files.o freeze.o group.o ielscf.o inform.o iounit.o keys.o mdstuf.o molcul.o moldyn.o mpole.o polar.o rgddyn.o rigid.o sizes.o stodyn.o units.o uprior.o usage.o
mdrest.o: atomid.o atoms.o bound.o group.o inform.o iounit.o mdstuf.o moldyn.o rgddyn.o sizes.o units.o
mdsave.o: atomid.o atoms.o bound.o boxes.o files.o group.o inform.o iounit.o mdstuf.o moldyn.o mpole.o output.o polar.o potent.o rgddyn.o sizes.o socket.o titles.o units.o
mdstat.o: atoms.o bath.o bound.o boxes.o inform.o inter.o iounit.o limits.o mdstuf.o molcul.o sizes.o units.o usage.o warp.o
mdstuf.o:
mechanic.o: ctran.o inform.o iounit.o limits.o potent.o vdwpot.o
merck.o: sizes.o
merge.o: atomid.o atoms.o couple.o iounit.o refer.o sizes.o
minima.o:
minimize.o: atoms.o files.o inform.o iounit.o keys.o scales.o sizes.o usage.o
minirot.o: files.o inform.o iounit.o keys.o math.o omega.o scales.o sizes.o zcoord.o
minrigid.o: files.o group.o inform.o iounit.o keys.o math.o output.o rigid.o sizes.o
molcul.o:
moldyn.o:
molecule.o: atomid.o atoms.o couple.o molcul.o sizes.o
molxyz.o: files.o iounit.o titles.o
moment.o:
moments.o: atomid.o atoms.o bound.o charge.o dipole.o limits.o moment.o mpole.o polar.o potent.o rigid.o sizes.o solute.o units.o usage.o
monte.o: atoms.o files.o inform.o iounit.o omega.o output.o sizes.o units.o usage.o zcoord.o
mplpot.o:
mpole.o:
mrecip.o:
mutant.o:
mutate.o: atomid.o atoms.o charge.o inform.o iounit.o katoms.o keys.o mpole.o mutant.o polar.o potent.o sizes.o
nblist.o: atoms.o bound.o boxes.o cell.o charge.o ctran.o iounit.o light.o limits.o mpole.o neigh.o potent.o sizes.o vdw.o
neigh.o:
newton.o: atoms.o files.o inform.o iounit.o keys.o sizes.o usage.o
newtrot.o: files.o hescut.o inform.o iounit.o keys.o math.o omega.o sizes.o zcoord.o
nextarg.o: argue.o
nexttext.o:
nonpol.o:
nose.o: atomid.o atoms.o bath.o boxes.o freeze.o mdstuf.o moldyn.o sizes.o units.o usage.o virial.o
nspline.o:
nucleic.o: atoms.o couple.o files.o group.o inform.o iounit.o katoms.o math.o molcul.o nucleo.o output.o potent.o resdue.o restrn.o rigid.o sequen.o sizes.o titles.o usage.o
nucleo.o: sizes.o
number.o: inform.o iounit.o
numeral.o:
numgrad.o: atoms.o sizes.o
ocvm.o: inform.o iounit.o keys.o linmin.o math.o minima.o output.o potent.o scales.o
omega.o:
opbend.o:
opdist.o:
openend.o:
openmp.o:
optimize.o: atoms.o files.o inform.o iounit.o keys.o scales.o sizes.o usage.o
optirot.o: files.o inform.o iounit.o keys.o math.o omega.o scales.o sizes.o zcoord.o
optrigid.o: files.o group.o inform.o iounit.o keys.o math.o output.o rigid.o
optsave.o: atoms.o files.o iounit.o math.o omega.o output.o scales.o sizes.o socket.o usage.o zcoord.o
orbital.o: atomid.o atoms.o bndstr.o couple.o iounit.o keys.o piorbs.o potent.o sizes.o tors.o
orbits.o:
orient.o: atomid.o atoms.o group.o math.o rigid.o sizes.o
orthog.o:
output.o:
overlap.o: units.o
params.o:
path.o: align.o atomid.o atoms.o files.o inform.o iounit.o linmin.o minima.o output.o paths.o sizes.o
paths.o:
pbstuf.o:
pdb.o:
pdbxyz.o: atomid.o atoms.o couple.o fields.o files.o inform.o iounit.o katoms.o pdb.o resdue.o sequen.o sizes.o titles.o
phipsi.o: sizes.o
picalc.o: atomid.o atoms.o bndstr.o couple.o inform.o iounit.o orbits.o piorbs.o pistuf.o sizes.o tors.o units.o
piorbs.o:
pistuf.o:
pitors.o:
pme.o:
pmestuff.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
pmpb.o: iounit.o
polar.o:
polarize.o: atoms.o inform.o iounit.o molcul.o mpole.o polar.o polgrp.o polpot.o potent.o sizes.o units.o
poledit.o: atomid.o atoms.o couple.o files.o iounit.o keys.o kpolr.o mpole.o polar.o polgrp.o polpot.o potent.o sizes.o units.o
polgrp.o:
polpot.o:
polymer.o: atoms.o bndstr.o bound.o boxes.o iounit.o keys.o sizes.o
potent.o:
potential.o: atomid.o atoms.o charge.o chgpot.o dipole.o files.o inform.o iounit.o katoms.o keys.o math.o minima.o moment.o mpole.o neigh.o output.o polar.o potent.o potfit.o ptable.o refer.o sizes.o titles.o units.o
potfit.o: sizes.o
precis.o:
precise.o:
pressure.o: atomid.o atoms.o bath.o bound.o boxes.o group.o iounit.o math.o mdstuf.o molcul.o moldyn.o sizes.o units.o usage.o virial.o
prmedit.o: angpot.o bndpot.o iounit.o math.o params.o sizes.o urypot.o vdwpot.o
prmkey.o: angpot.o bndpot.o chgpot.o ctran.o fields.o mplpot.o polpot.o potent.o rxnpot.o torpot.o urypot.o vdwpot.o
promo.o: iounit.o
protein.o: atomid.o atoms.o couple.o files.o group.o inform.o iounit.o katoms.o math.o molcul.o output.o phipsi.o potent.o resdue.o restrn.o rigid.o sequen.o sizes.o titles.o usage.o
prtdyn.o: atoms.o boxes.o files.o group.o mdstuf.o moldyn.o rgddyn.o sizes.o titles.o
prterr.o: files.o output.o
prtint.o: atomid.o atoms.o files.o inform.o sizes.o titles.o zclose.o zcoord.o
prtmol2.o: atomid.o atoms.o bndstr.o couple.o files.o iounit.o sizes.o titles.o
prtpdb.o: files.o pdb.o sequen.o sizes.o titles.o
prtprm.o: angpot.o bndpot.o chgpot.o ctran.o fields.o kanang.o kangs.o kantor.o katoms.o kbonds.o kchrge.o kdipol.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o mplpot.o polpot.o sizes.o urypot.o vdwpot.o
prtseq.o: files.o sequen.o sizes.o
prtxyz.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o sizes.o titles.o
pss.o: atoms.o files.o hescut.o inform.o iounit.o math.o omega.o refer.o sizes.o tree.o warp.o zcoord.o
pssrigid.o: atoms.o files.o group.o inform.o iounit.o math.o minima.o molcul.o refer.o rigid.o sizes.o warp.o
pssrot.o: atoms.o files.o inform.o iounit.o math.o minima.o omega.o refer.o sizes.o warp.o zcoord.o
ptable.o:
qmstuf.o:
qrfact.o:
quatfit.o: align.o sizes.o
radial.o: argue.o atomid.o atoms.o bound.o boxes.o files.o inform.o iounit.o limits.o math.o molcul.o potent.o sizes.o
random.o: inform.o iounit.o keys.o math.o
rattle.o: atomid.o atoms.o freeze.o group.o inform.o iounit.o moldyn.o sizes.o units.o usage.o virial.o
readdyn.o: atoms.o boxes.o files.o group.o iounit.o mdstuf.o moldyn.o rgddyn.o sizes.o
readgau.o: ascii.o iounit.o qmstuf.o sizes.o units.o
readgdma.o: atomid.o atoms.o dma.o files.o iounit.o mpole.o sizes.o units.o
readint.o: atomid.o atoms.o files.o inform.o iounit.o sizes.o titles.o zclose.o zcoord.o
readmol.o: atomid.o atoms.o couple.o files.o iounit.o ptable.o sizes.o titles.o
readmol2.o: atomid.o atoms.o couple.o files.o iounit.o ptable.o sizes.o titles.o
readpdb.o: files.o inform.o iounit.o pdb.o resdue.o sequen.o titles.o
readprm.o: ctran.o fields.o iounit.o kanang.o kangs.o kantor.o katoms.o kbonds.o kchrge.o kdipol.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o merck.o params.o sizes.o
readseq.o: files.o iounit.o resdue.o sequen.o
readxyz.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o iounit.o sizes.o titles.o
refer.o: sizes.o
replica.o: bound.o boxes.o cell.o inform.o iounit.o
resdue.o:
respa.o: atomid.o atoms.o freeze.o ielscf.o limits.o moldyn.o polar.o potent.o sizes.o units.o usage.o virial.o
restrn.o:
rgddyn.o:
rgdstep.o: atomid.o atoms.o bound.o group.o iounit.o rgddyn.o sizes.o units.o virial.o
rigid.o:
ring.o:
rings.o: angbnd.o atoms.o bitor.o bndstr.o couple.o inform.o iounit.o ring.o tors.o
rmsfit.o: align.o sizes.o
rotbnd.o:
rotlist.o: atoms.o couple.o iounit.o molcul.o rotbnd.o sizes.o zclose.o
rotpole.o: atoms.o mpole.o sizes.o
rxnfld.o:
rxnpot.o:
saddle.o: atoms.o inform.o iounit.o keys.o linmin.o minima.o sizes.o syntrn.o titles.o zcoord.o
scales.o: sizes.o
scan.o: atoms.o files.o inform.o iounit.o math.o minima.o omega.o output.o sizes.o zcoord.o
sdstep.o: atomid.o atoms.o bath.o couple.o freeze.o kvdws.o math.o mdstuf.o moldyn.o sizes.o stodyn.o units.o usage.o virial.o
search.o: linmin.o math.o sizes.o
sequen.o: sizes.o
server.o:
shakeup.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o bound.o couple.o freeze.o keys.o math.o molcul.o ring.o sizes.o usage.o
shunt.o:
sigmoid.o:
simplex.o: iounit.o keys.o minima.o
sizes.o:
sktstuff.o: atomid.o atoms.o charge.o couple.o deriv.o fields.o files.o inform.o iounit.o keys.o moldyn.o mpole.o polar.o potent.o sizes.o socket.o
sniffer.o: atoms.o files.o inform.o iounit.o linmin.o math.o minima.o output.o scales.o sizes.o usage.o
socket.o:
solute.o:
sort.o:
spacefill.o: atomid.o atoms.o files.o inform.o iounit.o kvdws.o math.o sizes.o usage.o
spectrum.o: files.o iounit.o math.o units.o
square.o: inform.o iounit.o keys.o minima.o sizes.o
stodyn.o:
strbnd.o:
strtor.o:
suffix.o: ascii.o
superpose.o: align.o atomid.o atoms.o bound.o files.o inform.o iounit.o sizes.o titles.o
surface.o: atoms.o inform.o iounit.o math.o sizes.o usage.o
surfatom.o: atoms.o iounit.o math.o sizes.o
switch.o: limits.o nonpol.o shunt.o
sybylxyz.o: files.o iounit.o titles.o
syntrn.o:
tarray.o:
temper.o: atomid.o atoms.o bath.o group.o ielscf.o mdstuf.o molcul.o moldyn.o rgddyn.o sizes.o units.o usage.o
testgrad.o: atoms.o deriv.o energi.o files.o inform.o inter.o iounit.o sizes.o solute.o usage.o
testhess.o: atoms.o files.o hescut.o inform.o iounit.o sizes.o usage.o
testpair.o: atoms.o deriv.o energi.o inform.o iounit.o light.o limits.o neigh.o polpot.o potent.o sizes.o tarray.o vdwpot.o
testpol.o: atoms.o bound.o inform.o iounit.o limits.o minima.o polar.o polpot.o potent.o rigid.o sizes.o units.o usage.o
testrot.o: domega.o energi.o inform.o iounit.o math.o omega.o sizes.o zcoord.o
timer.o: atoms.o hescut.o inform.o iounit.o limits.o sizes.o
timerot.o: iounit.o limits.o omega.o sizes.o
titles.o:
tncg.o: atoms.o hescut.o inform.o iounit.o keys.o linmin.o math.o minima.o output.o piorbs.o potent.o sizes.o
torphase.o:
torpot.o:
torque.o: atoms.o deriv.o mpole.o sizes.o
tors.o:
torsfit.o: atomid.o atoms.o files.o inform.o iounit.o keys.o ktorsn.o math.o output.o potent.o qmstuf.o restrn.o scales.o sizes.o tors.o usage.o
torsions.o: atoms.o bndstr.o couple.o iounit.o sizes.o tors.o
tortor.o:
tree.o:
trimtext.o:
unitcell.o: bound.o boxes.o iounit.o keys.o
units.o:
uprior.o:
urey.o:
urypot.o:
usage.o:
usolve.o:
valence.o: angbnd.o angpot.o atomid.o atoms.o bndpot.o bndstr.o couple.o files.o hescut.o inform.o iounit.o kangs.o kbonds.o keys.o kopbnd.o kstbnd.o ktorsn.o kurybr.o kvdws.o linmin.o math.o minima.o opbend.o output.o potent.o qmstuf.o scales.o sizes.o strbnd.o torpot.o tors.o units.o urey.o urypot.o usage.o valfit.o vdwpot.o
valfit.o:
vdw.o:
vdwpot.o:
verlet.o: atomid.o atoms.o freeze.o ielscf.o moldyn.o polar.o sizes.o units.o usage.o
version.o: iounit.o output.o
vibbig.o: atomid.o atoms.o bound.o couple.o files.o hescut.o hessn.o inform.o iounit.o keys.o limits.o mpole.o potent.o rigid.o sizes.o units.o usage.o vdw.o vdwpot.o vibs.o
vibrate.o: atomid.o atoms.o files.o hescut.o iounit.o math.o sizes.o units.o usage.o
vibrot.o: iounit.o omega.o sizes.o
vibs.o:
virial.o:
volume.o: atoms.o iounit.o math.o sizes.o
warp.o:
xtalfit.o: atomid.o atoms.o bound.o boxes.o charge.o dipole.o energi.o files.o fracs.o inform.o iounit.o kvdws.o limits.o math.o molcul.o mpole.o polar.o potent.o sizes.o vdw.o vdwpot.o xtals.o
xtalmin.o: atoms.o boxes.o files.o inform.o iounit.o keys.o math.o scales.o sizes.o
xtals.o:
xyzatm.o: atoms.o inform.o iounit.o math.o sizes.o
xyzedit.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o iounit.o limits.o linmin.o math.o minima.o molcul.o neigh.o output.o potent.o ptable.o refer.o scales.o sizes.o titles.o units.o usage.o
xyzint.o: files.o iounit.o titles.o
xyzpdb.o: atomid.o atoms.o couple.o fields.o files.o inform.o molcul.o pdb.o resdue.o sequen.o sizes.o
xyzsybyl.o: files.o iounit.o sizes.o titles.o
zatom.o: angbnd.o atomid.o atoms.o bndstr.o fields.o iounit.o kangs.o katoms.o kbonds.o sizes.o zclose.o zcoord.o
zclose.o: sizes.o
zcoord.o: sizes.o
