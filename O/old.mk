#   OPTIM: A program for optimizing geometries and calculating reaction pathways
#   Copyright (C) 1999-2006 David J. Wales
#   This file is part of OPTIM.
#
#   OPTIM is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   OPTIM is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
OBJS = potential.o key.o commons.o modtwoend.o modhess.o modneb.o modamber.o modamber2.o modamber9.o grouprotation.o\
modunres.o modguess.o modmec.o syminf.o vecck.o zwk.o \
getparams.o adm.o amb_natinterns.o capsid.o axdiff.o axpairs.o aziz.o c60diff.o c60p.o beig.o bfgsts.o \
dimer.o guesspath.o unguesspath.o \
charmmBildc.o chiralhyd.o dtrap.o dumpit.o dumpp.o eig.o eigensort.o emie.o escp.o etrap.o fetchz.o geopt.o gmetry.o hessout.o \
ions.o JM.o keywords.o h2o.o latmin.o ljdiff.o ljpdiff.o mdiff.o mied.o miel.o mlatmin.o \
morse.o rotd.o mpdiff.o msdiff.o mslatmin.o oldneb.o pertable.o shifth.o dprand.o \
rotcon.o rotm.o scdiff.o scl.o secdiag.o siaz.o sortc.o sortxyz.o symmetry.o tcheck.o vdump.o \
efol.o utils.o ljpbin.o ljpshift.o SW.o ctest.o double.o Clatmin.o input.o cpmdlatmin.o ljpkob.o \
changep.o twoend.o dzugutov.o Zetterling.o rad.o odesd.o g2special.o PV.o PachecoC60.o xmylbfgs.o \
mylbfgs.o 2Dfunc.o OPTIM.o path.o connect.o mindist.o repelsp.o TIPmodes.o caardiff.o mycpu_time.o \
g46merdiff.o p46merdiff.o MB.o c10.o ectrap.o dctrap.o welch.o tosifumi.o tight.o dftb.o \
fd.o ljms.o meccano.o unmeccano.o porfuncs.o modcharmm.o modmxatms.o rigidfuncs.o tip.o \
msevb_common.o msevb_interfaces.o checkedtrig.o msevb.o fmsevb.o drvmsevb.o cleanmemory.o Natb.o \
header.o wc.o binaryio.o BLN.o morph.o newmindist.o thomson.o minpermdist.o minperm.o inertia.o \
greatcirc.o GoOptim.o cubspl.o specialfuncts.o gsdata.o cubsplstring.o dqag.o \
intcommons.o intcoords.o SiO2.o stock.o bhinterp.o Z2faster.o p4diff.o Gupta.o myorient.o hybridmin.o \
bbrsdm.o projh.o bisect.o dbpg.o dbptd.o lwotpbox.o lwotpgh.o gb.o gbbox.o gbd.o rbperm.o rigidb.o msstock.o \
newcapsid.o newtip.o pahagh.o pyg.o rbinertia.o stockghaa.o gay-berne.o growstring.o ambnatintern_extras.o minpermrb.o \
amhglobals.o Ackland_wrapper.o Ackland_metals.o DF1.o bulkmindist.o SD.o finite_differences.o minpermdistrbcom.o ttm3f.o \
make_conpot.o qtip4pf.o qspcfw.o intlbfgs.o congrad.o intlbfgslj.o intgradlj.o H2O.pjt2.o VRTMCY5f.o \
EYtrap.o patchyd.o asaoos.o

# Set up preprocessing for .F files
 CPP = /lib/cpp
 CPFLAGS = -traditional -P
# 

GENF90FILES = porfuncs.f90 header.f90

CHDUM = chdummy.o
UNDUM = unresdummy.o
AMDUM = amberdummy.o
AMHDUM = amhdummy.o
JBDUM = libjbdummy.a

VPATH = NEB:CONNECT:CHARMM:AMH:.
LDFLAGS = -L.

#### libraries
LIBS = -lnn  -lnc libmyblas.a libmylapack.a  libmyblas.a  libmyblas.a -lnc
# LIBS = -lmylapack -lmyblas -L/usr/local/lib
# LIBS = -lnn -lnc -lATLASf77blas -lATLAScblas -lATLASatlas -lATLASlapack

#### sf344> AMBER 9 / NAB interface
AMB9DUM = amber9dummy.o
AMB9SRC = ../../AMBER/src/sander
NABSRC = ../../NAB/src
#### end AMBER 9

#
#### AMBER
EXTRASAMBER = amber.o
# EXTRASAMH = amhglobals.o
EXTRASAMH = 
#### end AMBER
#

ifeq (${CTYPE},C31)
##### CHARMM 31
EXTRAS = myblas.o mylapack.o
CHOBJS31 = $(LIBDIR31)/help.o $(LIBDIR31)/iniall.o $(LIBDIR31)/miscom.o $(LIBDIR31)/usersb.o
CHLIBS31 = $(LIBDIR31)/adumb.a \
$(LIBDIR31)/flucq.a $(LIBDIR31)/cadint.a $(LIBDIR31)/cheq.a $(LIBDIR31)/cff.a $(LIBDIR31)/correl.a $(LIBDIR31)/dimb.a \
$(LIBDIR31)/emap.a $(LIBDIR31)/dynamc.a $(LIBDIR31)/energy.a $(LIBDIR31)/gamint.a $(LIBDIR31)/gukint.a \
$(LIBDIR31)/gener.a $(LIBDIR31)/image.a $(LIBDIR31)/io.a $(LIBDIR31)/machdep.a $(LIBDIR31)/manip.a $(LIBDIR31)/mbond.a \
$(LIBDIR31)/mc.a $(LIBDIR31)/minmiz.a $(LIBDIR31)/misc.a $(LIBDIR31)/mmff.a $(LIBDIR31)/molvib.a $(LIBDIR31)/nbonds.a \
$(LIBDIR31)/pert.a $(LIBDIR31)/quantum.a $(LIBDIR31)/rxncor.a $(LIBDIR31)/shapes.a $(LIBDIR31)/solvation.a \
$(LIBDIR31)/util.a $(LIBDIR31)/vibran.a libcharmm.a 
##### end CHARMM 31

else

##### CHARMM 35
EXTRAS = myblas.o mylapack.o
CHOBJS35 = $(LIBDIR35)/help.o $(LIBDIR35)/iniall.o $(LIBDIR35)/miscom.o $(LIBDIR35)/usersb.o
CHLIBS35 = $(LIBDIR35)/adumb.a \
$(LIBDIR35)/flucq.a $(LIBDIR35)/cadint.a $(LIBDIR35)/cheq.a $(LIBDIR35)/cff.a $(LIBDIR35)/correl.a $(LIBDIR35)/dimb.a \
$(LIBDIR35)/emap.a $(LIBDIR35)/dynamc.a $(LIBDIR35)/energy.a $(LIBDIR35)/gamint.a $(LIBDIR35)/gukint.a \
$(LIBDIR35)/gener.a $(LIBDIR35)/image.a $(LIBDIR35)/io.a $(LIBDIR35)/machdep.a $(LIBDIR35)/manip.a $(LIBDIR35)/mbond.a \
$(LIBDIR35)/mc.a $(LIBDIR35)/minmiz.a $(LIBDIR35)/misc.a $(LIBDIR35)/mmff.a $(LIBDIR35)/molvib.a $(LIBDIR35)/nbonds.a \
$(LIBDIR35)/pert.a $(LIBDIR35)/quantum.a $(LIBDIR35)/rxncor.a $(LIBDIR35)/shapes.a $(LIBDIR35)/solvation.a \
$(LIBDIR35)/util.a $(LIBDIR35)/vibran.a libcharmm.a
##### end CHARMM 35
endif


#
##### UNRES
SRC = /export/home/wales/unres/src
EXTRASUNRES=$(SRC)/unres_subroutine.o $(SRC)/arcos.o $(SRC)/cartprint.o $(SRC)/chainbuild.o \
$(SRC)/convert.o $(SRC)/initialize_p.o $(SRC)/matmult.o $(SRC)/readrtns_p.o \
$(SRC)/parmread.o $(SRC)/gen_rand_conf.o $(SRC)/printmat.o $(SRC)/map.o \
$(SRC)/pinorm.o $(SRC)/randgens.o $(SRC)/rescode.o $(SRC)/intcor.o $(SRC)/timing.o \
$(SRC)/misc.o $(SRC)/cartder.o $(SRC)/checkder_p.o $(SRC)/energy_p.o $(SRC)/gradient_p.o \
$(SRC)/minimize_p.o $(SRC)/sumsld.o $(SRC)/cored.o $(SRC)/rmdd.o $(SRC)/geomout.o \
$(SRC)/readpdb.o $(SRC)/regularize.o $(SRC)/thread.o $(SRC)/fitsq.o $(SRC)/mcm.o \
$(SRC)/mc.o $(SRC)/bond_move.o $(SRC)/refsys.o $(SRC)/check_sc_distr.o \
$(SRC)/contact.o $(SRC)/djacob.o $(SRC)/entmcm.o $(SRC)/compare_s1.o \
intbfgsts.o intbeig.o intxmylbfgs.o intsecdiag.o \
unresoptim.o unresdump.o unrestransform.o unressetdihe.o unrestwist.o unresconnectsections.o \
unresguessts.o unrescalcdihe.o unresnebguessts.o 
##### end UNRES

###############################################################################################
# Add the OPTIM version number to OPTIM.f so it is printed in the output
DEFS=
DEFS+=-DSVNVERSION="`./version.sh`"
###################################### RULES AND TARGETS ######################################
.SUFFIXES: 
.SUFFIXES: .o .f .F .f90

.F.f:
	$(CPP) $(CPFLAGS) $(DEFS) $< > $@

.f90.o:
	$(FC) $(FFLAGS) ${SEARCH_PATH} -c $<
.f.o:
	$(FC) $(FFLAGS) ${SEARCH_PATH} -c $<

default:
	./build.csh help
${optim}: $(AMB9DUM) $(AMHDUM) $(JBDUM) $(OBJS) $(EXTRAS) $(CHDUM) $(UNDUM) $(AMDUM)
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o ${optim} $(EXTRAS) $(OBJS) $(CHDUM) $(UNDUM)  $(AMHDUM) $(AMDUM) $(AMB9DUM) $(JBDUM) \
	$(LDFLAGS) $(LIBS) $(LIBS)
	rm -f libmyblas.a    libmylapack.a

ifeq (${CTYPE},C31)
${coptim}: $(AMB9DUM) $(AMHDUM) $(JBDUM) $(OBJS) $(EXTRAS) $(UNDUM) $(AMDUM) libcharmm.a 
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o ${coptim} ${CHOBJS31} $(EXTRAS) $(UNDUM) $(AMDUM) $(AMHDUM) $(AMB9DUM) $(OBJS) $(JBDUM) \
	$(LDFLAGS) $(LIBS) $(CHLIBS31) ${CHLIBS31} ${CHLIBS31} 
	rm libmyblas.a    libmylapack.a

else

${coptim}: $(AMB9DUM) $(AMHDUM) $(JBDUM) $(OBJS) $(EXTRAS) $(UNDUM) $(AMDUM) libcharmm.a 
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o ${coptim} ${CHOBJS35} $(EXTRAS) $(UNDUM) $(AMDUM) $(AMHDUM) $(AMB9DUM) $(OBJS) $(JBDUM) \
	$(LDFLAGS) $(LIBS) $(CHLIBS35) ${CHLIBS35} ${CHLIBS35} 
	rm libmyblas.a    libmylapack.a
endif

${unoptim}: $(AMB9DUM) $(AMHDUM) $(JBDUM) $(OBJS) $(EXTRAS) $(EXTRASUNRES) $(CHDUM) $(AMDUM)
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o ${unoptim} $(EXTRAS) $(EXTRASUNRES) $(CHDUM) $(AMDUM) $(AMHDUM) $(AMB9DUM) $(OBJS) $(JBDUM) \
	$(LDFLAGS) $(LIBS)
	rm libmyblas.a    libmylapack.a

${amoptim}: $(AMB9DUM) $(AMHDUM) $(JBDUM) $(OBJS) $(EXTRAS) $(EXTRASAMBER) $(CHDUM) $(UNDUM) 
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o ${amoptim} $(EXTRAS) $(EXTRASAMBER) $(CHDUM) $(UNDUM) $(AMHDUM) $(AMB9DUM) $(OBJS) $(JBDUM) \
	$(LDFLAGS) $(LIBS)
	rm libmyblas.a    libmylapack.a

${amb9optim}: $(AMHDUM) $(EXTRAS) $(JBDUM) $(OBJS) $(CHDUM) $(AMDUM) $(UNDUM) libamber.a libnab.a 
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o ${amb9optim} $(CHDUM) $(UNDUM) $(AMHDUM) $(AMDUM) $(OBJS) $(JBDUM) $(EXTRAS) \
	$(LDFLAGS) $(LIBS) libamber.a libnab.a $(LIBS)
#	rm -f libmyblas.a    libmylapack.a modamber9.o modamber9.mod

${amhoptim}: $(AMB9DUM) $(JBDUM) $(OBJS) $(EXTRAS) $(AMHEXTRAS) $(CHDUM) $(AMDUM) $(UNDUM) libamh.a 
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o ${amhoptim} $(EXTRAS)  $(OBJS) $(AMHEXTRAS) $(CHDUM) $(UNDUM) $(AMDUM) $(AMB9DUM) $(JBDUM) \
	$(LDFLAGS) $(LIBS)  $(LIBS) libamh.a
	rm libmyblas.a    libmylapack.a 

${jboptim}: $(AMB9DUM) $(AMHDUM) libbowman.a $(OBJS) $(EXTRAS) $(CHDUM) $(UNDUM) $(AMDUM)
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o ${jboptim} $(EXTRAS) $(OBJS) $(CHDUM) $(UNDUM)  $(AMHDUM) $(AMDUM) $(AMB9DUM) \
	$(LDFLAGS) $(LIBS) libbowman.a Bowman/libs/libpes.a Bowman/libs/libpesd.a
	rm -f libmyblas.a    libmylapack.a

timing:
	rm -f $(OPROG) 
	$(FC) $(FFLAGS) -p $(LDFLAGS) -o $(OPROG) $(OBJS) $(EXTRAS) $(LIBS)
	rm libmyblas.a    libmylapack.a

clean:
	rm -f $(OBJS) *.mod $(EXTRAS) *.lst *.o pref.dat prefx.msg FOR021.DAT ${GENF90FILES} *.log lastargs lib*
	rm -f libbowman.a libjbdummy.a
	if test -d ../../BLAS ;  then cd ../../BLAS ; make clean ; fi
	if test -d ../../LAPACK ;  then cd ../../LAPACK ; make clean ; fi
	if test -d NEB ;  then cd NEB ; make clean ; fi
	if test -d CONNECT ;  then cd CONNECT ; make clean ; fi
	if test -d CHARMM ;  then cd CHARMM ; make clean ; fi
	if test -d CHARMM35 ;  then cd CHARMM35 ; make clean ; fi
	if test -d AMH ;  then cd AMH ; make clean ; fi
	if test -d $(AMB9SRC) ; then cd $(AMB9SRC) ; make clean ; fi
	if test -d $(NABSRC) ; then cd $(NABSRC) ; make clean ; fi
	if test -d Bowman/fake ; then cd Bowman/fake/ ; make clean ; fi
	if test -d Bowman ; then cd Bowman ; make clean ; fi
cleanexe:
	rm -r ../bin/*

cleanall:
	make clean
	make cleanexe
cleanamber:
	rm -f libmyblas.a libmylapack.a modamber9.o modamber9.mod
porfuncs.f90: porfuncs.csh
	./porfuncs.csh ${SWITCH} > porfuncs.f90
header.f90: header.csh
	./header.csh > header.f90

# SAT: no agressive optimizations for selected files to cut down on compile time
keywords.o: keywords.f
	${FC} ${NOOPT} ${SEARCH_PATH} -c keywords.f
libamber.a:
	cd $(AMB9SRC); make FC="${FC}" FFLAGS="${FFLAGS}" FREEFORMAT_FLAG="${FREEFORMAT_FLAG}" EXTRA_FLAGS="${EXTRA_FLAGS}" SRCDIR="$(CURDIR)"
libnab.a:
	cd $(NABSRC)/../ucpp-1.3/; make install; cd $(NABSRC); make install NABHOME="$(CURDIR)/../../NAB" LIBDIR="$(CURDIR)/../../NAB/lib";cp $(CURDIR)/../../NAB/lib/libnab.a $(CURDIR)/libnab.a
libnn.a: SAT-Ghost
	cd NEB; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}"
libamh.a: SAT-Ghost
	cd AMH; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}" 
libnc.a: SAT-Ghost
	cd CONNECT; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}"
libmyblas.a: SAT-Ghost
	cd ../../BLAS; make double FC="${FC}" FFLAGS="${FFLAGS}" BLAS_EXCLUDE_LIST="${BLAS_EXCLUDE_LIST}";\
	cp libmyblas.a ../OPTIM/source
libmylapack.a: SAT-Ghost
	cd ../../LAPACK; make selection FC="${FC}" FFLAGS="${FFLAGS}" NOOPT="${NOOPT}";\
	cp libmylapack.a ../OPTIM/source

ifeq (${CTYPE},C31)
libcharmm.a: SAT-Ghost
	cd CHARMM; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}" PREFLX="${PREFLX}" PREFDIR="${PREFDIR}" \
	CTYPE="${CTYPE}" FCMDIR=${FCMDIR} C31SRC="${C31SRC}" SRC31="${SRC31}"

else

libcharmm.a: SAT-Ghost
	cd CHARMM35; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}" PREFLX="${PREFLX}" PREFDIR="${PREFDIR}" \
	CTYPE="${CTYPE}" FCMDIR=${FCMDIR} C35SRC="${C35SRC}" SRC35="${SRC35}"

endif

libbowman.a:
	cd Bowman; make
libjbdummy.a:
	cd Bowman/fake; make FC="${FC}" FFLAGS="${FFLAGS}"
SAT-Ghost:

###################################### DEPENDENCIES ######################################
${optim} ${coptim} ${unoptim} ${amoptim} ${amb9optim} ${amhoptim} ${jboptim}: libnn.a libnc.a libmyblas.a libmylapack.a
${amb9optim}: libamber.a modamber9.o
# libamh.a: altpot_interfaces.o  amh_interfaces.o  globals_alt.o amhglobals.o
libamh.a: amhglobals.o
libnn.a: modcharmm.o key.o commons.o modunres.o porfuncs.o modmec.o modguess.o efol.o growstring.o gsdata.o intcommons.o intcoords.o amhglobals.o 
libnc.a: libnn.a key.o syminf.o modhess.o amhglobals.o
libamber.a: commons.o modamber9.o grouprotation.o
libcharmm.a: key.o commons.o modtwoend.o modneb.o modhess.o modcharmm.o modmxatms.o intcommons.o

2Dfunc.o:      commons.o modhess.o
Ackland_wrapper.o:      key.o
BLN.o:      key.o modhess.o
Clatmin.o:      commons.o
DF1.o:      key.o modhess.o
GoOptim.o:      key.o
JM.o:      modhess.o porfuncs.o
MB.o:      modhess.o commons.o
PV.o:      commons.o key.o
PachecoC60.o:      modhess.o
SW.o:      key.o modhess.o
SiO2.o:      porfuncs.o modhess.o
TIPmodes.o:      modhess.o
Zetterling.o:      porfuncs.o modhess.o key.o
adm.o:      porfuncs.o commons.o key.o
amber.o:      commons.o modamber.o modamber2.o modhess.o key.o
axdiff.o:      modhess.o
aziz.o:      modhess.o
beig.o:      commons.o key.o modneb.o modtwoend.o porfuncs.o
bfgsts.f:      commons.o key.o vecck.o zwk.o modcharmm.o modhess.o porfuncs.o modtwoend.o 
c10.o:        modhess.o
c60diff.o:      modhess.o
c60p.o:      modhess.o
caardiff.o:      modhess.o
capsid.o:      modhess.o
changep.o:      commons.o key.o modtwoend.o modamber.o modneb.o porfuncs.o
charmmBildc.o:      modamber9.o commons.o key.o
connect.o:      commons.o key.o modtwoend.o modneb.o modcharmm.o modunres.o porfuncs.o porfuncs.o syminf.o modhess.o modcharmm.o modunres.o porfuncs.o porfuncs.o porfuncs.o modneb.o key_neb.o libnn.a modcharmm.o modunres.o efol.o porfuncs.o modcharmm.o modamber9.o 

cpmdlatmin.o:      key.o porfuncs.o
ctest.o:      modhess.o
dctrap.o:      key.o modhess.o
dftb.o:      key.o modhess.o
diis.o:      key.o gdiis.o
dimer.o:       modhess.o commons.o modhess.o commons.o
double.o:      modhess.o
dtrap.o:      modhess.o
dumpit.o:      commons.o key.o
dumpp.o:      commons.o key.o  modcharmm.o porfuncs.o
dzugutov.o:      modhess.o
ectrap.o:      key.o
escp.o:      modhess.o
fd.o:      modhess.o key.o
fetchz.o:      commons.o key.o modtwoend.o syminf.o modneb.o modcharmm.o modamber9.o modunres.o libnc.a modguess.o modmec.o porfuncs.o msevb_common.o binaryio.o amhglobals.o
g2special.o:      modhess.o key.o porfuncs.o
g46merdiff.o:        modhess.o
geopt.o:      commons.o key.o modtwoend.o modhess.o modcharmm.o porfuncs.o binaryio.o efol.o modneb.o amhglobals.o modhess.o modhess.o modhess.o modcharmm.o 
getparams.o:      commons.o porfuncs.o header.o key.o modamber9.o
gmetry.o:      key.o commons.o
greatcirc.o:      commons.o key.o modtwoend.o modhess.o zwk.o modunres.o modcharmm.o porfuncs.o
h2o.o:      modhess.o
hessout.o:      modhess.o
hybridmin.o:      commons.o key.o vecck.o zwk.o modcharmm.o modhess.o porfuncs.o
inertia.o:      commons.o key.o 
input.o:      porfuncs.o
intbeig.o:      commons.o key.o modneb.o modtwoend.o
intbfgsts.o:      commons.o key.o vecck.o zwk.o modunres.o modhess.o
bfgsts.o:      commons.o key.o  vecck.o zwk.o modcharmm.o modhess.o modtwoend.o porfuncs.o
intsecdiag.o:      commons.o key.o
intxmylbfgs.o:      commons.o key.o modtwoend.o
ions.o:      modhess.o
keywords.o:      commons.o key.o modmec.o modtwoend.o modamber.o modamber9.o modneb.o modcharmm.o modunres.o key_neb.o libnc.a modmec.o modguess.o porfuncs.o gay-berne.o msevb_common.o wc.o binaryio.o gsdata.o cubsplstring.o intcommons.o amhglobals.o modmxatms.o
ljdiff.o:      modhess.o
ljms.o:      modhess.o
ljpbin.o:      key.o modhess.o
ljpdiff.o:      key.o modhess.o
ljpkob.o:      key.o modhess.o
ljpshift.o:      key.o modhess.o
mdiff.o:      modhess.o
mied.o:      modhess.o
mindist.o:      key.o porfuncs.o
morph.o:      commons.o key.o zwk.o modtwoend.o modcharmm.o porfuncs.o
morse.o:      modhess.o
mpdiff.o:      modhess.o key.o
msdiff.o:      modhess.o
mylbfgs.o:      commons.o key.o modtwoend.o modhess.o zwk.o modunres.o modcharmm.o modamber9.o porfuncs.o specialfuncts.o
myorient.o:      key.o
odesd.o:      commons.o key.o modhess.o porfuncs.o
oldneb.o:      commons.o modtwoend.o modneb.o key.o modhess.o porfuncs.o commons.o modcharmm.o
p46merdiff.o:        modhess.o
p4diff.o:      modhess.o
path.o:      commons.o key.o syminf.o modcharmm.o modunres.o modhess.o efol.o porfuncs.o amhglobals.o
pertable.o:      commons.o key.o
potential.o:      commons.o key.o modhess.o modcharmm.o porfuncs.o modamber9.o commons.o modhess.o
projh.o:      modhess.o
rad.o:      commons.o key.o
rbinertia.o: commons.o key.o
repelsp.o:      commons.o key.o
rotd.o:      modhess.o
scdiff.o:      modhess.o
secdiag.o:      commons.o key.o modcharmm.o porfuncs.o
shifth.o:      key.o modhess.o
sortxyz.o:      commons.o
stock.o:      key.o modhess.o commons.o vecck.o modhess.o
symmetry.o:      commons.o key.o syminf.o
tcheck.o:      commons.o
thomson.o:      modhess.o commons.o key.o 
tight.o:      key.o modhess.o
tip.o:      commons.o key.o modhess.o
tosifumi.o:      porfuncs.o modhess.o
twoend.o:      porfuncs.o commons.o key.o modtwoend.o zwk.o modcharmm.o modhess.o
unrescalcdihe.f:      commons.o modunres.o commons.o modunres.o modunres.o commons.o
unresconnectsections.o:      commons.o modunres.o commons.o key.o modtwoend.o modunres.o commons.o modtwoend.o modunres.o
unresdummy.o:      commons.o
unresdump.o:      commons.o
unresguessts.o:      commons.o key.o modtwoend.o modunres.o commons.o modtwoend.o modunres.o
unresnebguessts.o:      commons.o key.o modunres.o key_neb.o commons.o modunres.o
unresoptim.o:      commons.o modunres.o commons.o modunres.o commons.o modunres.o commons.o modunres.o
unressetdihe.o:      modunres.o
unrestransform.o:      modunres.o commons.o modhess.o modunres.o modunres.o
unrestwist.o:      commons.o modunres.o commons.o modunres.o commons.o modunres.o
utils.o:      commons.o
vdump.o:      key.o modhess.o porfuncs.o
welch.o:      porfuncs.o modhess.o
xmylbfgs.o:      commons.o key.o modtwoend.o porfuncs.o
amb_natinterns.o:      commons.o modamber9.o intcommons.o 
ambnatintern_extras.o:      modamber9.o intcommons.o modmxatms.o specialfuncts.o commons.o key.o 
bbrsdm.o:      commons.o key.o modcharmm.o
bhinterp.o: key.o modcharmm.o modamber9.o porfuncs.o commons.o efol.o specialfuncts.o
binaryio.o:          key.o commons.o porfuncs.o commons.o
bisect.o: key.o modcharmm.o modamber9.o porfuncs.o commons.o
checkedtrig.o:         msevb_common.o
chiralhyd.o:      commons.o modamber9.o intcommons.o commons.o modamber9.o key.o 
cleanmemory.o:  msevb_common.o
cubsplstring.o:  gsdata.o
dbpg.o:      modhess.o commons.o key.o commons.o 
dbptd.o:     modhess.o commons.o key.o commons.o 
drvmsevb.o:	commons.o msevb_common.o msevb_interfaces.o msevb_interfaces.o modhess.o
efol.o:      commons.o key.o syminf.o modhess.o zwk.o porfuncs.o modcharmm.o
fmsevb.o:	commons.o msevb_common.o msevb_interfaces.o
gay-berne.o: key.o commons.o commons.o modhess.o
gb.o:      modhess.o key.o
gbbox.o:      modhess.o commons.o key.o
gbd.o:      modhess.o key.o
growstring.o:  commons.o key.o cubsplstring.o gsdata.o specialfuncts.o intcommons.o intcoords.o
gsdata.o:  commons.o intcommons.o key.o
guesspath.o: key.o modcharmm.o modguess.o modunres.o porfuncs.o
intcommons.o:  modamber9.o
intcoords.o:  modcharmm.o commons.o intcommons.o key.o specialfuncts.o modamber9.o porfuncs.o
lwotpbox.o: modhess.o commons.o key.o 
lwotpgh.o: modhess.o commons.o key.o 
meccano.o: commons.o porfuncs.o key.o modneb.o modtwoend.o modmec.o key_neb.o commons.o modunres.o modcharmm.o 
meccano.o: key.o
minpermdist.o: key.o modcharmm.o modamber9.o intcommons.o intcoords.o commons.o
minpermrb.o: key.o modcharmm.o intcommons.o commons.o 
msevb.o:     commons.o msevb_common.o msevb_interfaces.o
msevb_interfaces.o:       msevb_common.o
msstock.o:     modhess.o commons.o key.o
mycpu_time.o: porfuncs.o key.o
newcapsid.o:      commons.o 
newmindist.o: commons.o key.o
newtip.o:      commons.o key.o modhess.o
pahagh.o: modhess.o commons.o key.o
patchyd.o: modhess.o commons.o
asaoos.o: modhess.o commons.o
pyg.o:    modhess.o commons.o key.o
rbperm.o: commons.o key.o
minpermdistrbcom.o: commons.o key.o
rigidb.o: key.o modhess.o commons.o
specialfuncts.o:    commons.o
stockghaa.o: modhess.o commons.o key.o 
unguesspath.o: key.o modcharmm.o modguess.o modunres.o porfuncs.o
unmeccano.o: commons.o key.o modneb.o modunres.o modmec.o key_neb.o commons.o modtwoend.o zwk.o modcharmm.o modneb.o modmec.o porfuncs.o modhess.o
OPTIM.f: OPTIM.F
OPTIM.o: commons.o key.o modtwoend.o modhess.o vecck.o zwk.o modcharmm.o modunres.o libnn.a libnc.a modguess.o modmec.o porfuncs.o charutils.o intcommons.o intcoords.o 
bulkmindist.o:  key.o
fetchz.o: SD.o 
#potential.o: SD.o finite_differences.o libbowman.a
potential.o: SD.o finite_differences.o ttm3f.o VRTMCY5f.o H2O.pjt2.o
make_conpot.o: key.o
grouprotation.o: commons.o
intlbfgs.o: key.o porfuncs.o commons.o intcommons.o amhglobals.o efol.o
intlbfgslj.o: key.o porfuncs.o commons.o intcommons.o efol.o
congrad.o: key.o commons.o 
intgradlj.o: key.o commons.o 
EYtrap.o: key.o
