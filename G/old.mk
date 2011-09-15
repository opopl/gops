# Introductory remarks {{{
# 
# Default target gives GMIN
# binaries are now produced in the bin directory!
# 'make chgmin' produces CGMIN
# 'make AMBGMIN' produces AMBGMIN
# 'make clean' removed all object and library files, and any binaries in the bin directory
#
#  DONT FORGET - mpirun
#              - add MPI keyword to data to run with mpirun
#              - change CHARMM library to the appropriate version
#
#  64 bit CGMIN openmpi executable works with BSPT LJ31. 
#  
#  Tried linking to CHARMM compiled with 64 bit openmpi. Reduced optimisation options.
#  BSPT metenk run works fine on one node and more than one node when linked to openmpi charmm build.
#  Also works when linked to 64 bit charmm library compiled with pgf90.
#  -C flag cannot be used with any CHARMM build due to CHARMM intentional out of bounds code.
#
#  pgf90 7.1.6 works for 64 bit lam and openmpi
#  pgf90 7.1.6 works for 64 bit openmpi with CHARMM library ~/svn/CHARMM35
#        for -Mextend -O3 -Munroll -Mnoframe flags in 
#        ~/svn/CHARMM35/build/UNX/Makefile_gnu and modules
#  pgi/64/7.1/6 and mpi/openmpi/64/pgi71/1.2.6
#
#  pgf90 7.1.6 64 bit lam cannot be tested due to a dependency
#
#  openmpi works with ifort on mek-quake:
#  ifort_em64t/10.0.025 mpi/openmpi/64/intel10/1.2.6 icc_em64t/10.0.023
#
#  openmpi fails on clust with module:
#
# ifort/64/10.1/015 icc/64/10.1/015 mpi/openmpi/64/intel10/1.2.6
#
# }}}
# Definitions {{{

GPROG =	../bin/GMIN
CPROG = ../bin/CGMIN
AMHPROG = ../bin/AMHGMIN
AMBPROG = ../bin/AMBGMIN

OBJS1 = porfuncs.o commons.o countatoms.o modamb.o modcharmm.o modmxatms.o modhess.o modamber9.o grouprotation.o \
	operations.o
OBJS2 =	centre.o finalio.o modconsts_trans_97.o modconsts.o dist.o\
	io1.o keyword.o main.o mc.o mcruns.o morse.o \
	potential.o quench.o rad.o dprand.o saveit.o seed.o \
	sort.o sort2.o sort3.o sort4.o takestep.o mycpu_time.o trans.o \
	finalq.o symmetry.o symmetrycsm.o ptgrp.o eigsrt.o SiSW.o taboo.o reseed.o newinertia.o supermc.o \
	tosifumi.o ortho.o compress.o mylbfgs.o input.o ddfpmin.o dlnsrch.o cgmin.o linmin.o \
	brent.o mnbrak.o dbrent.o f1dim.o zwischen.o hsmove.o PachecoC60.o AT.o EAMLJ_sub.o \
	Pbglue.o wenzel.o odesd.o capsid.o rigidfuncs.o tip.o pah.o strand.o \
	SW.o qmod.o ljpbin.o fdm.o dftb.o ljpshift.o dzugutov.o ljcoulomb.o \
	fd.o fedor.o welch.o BLJcluster.o stock.o Farkas.o getorbits.o \
	sc.o Zetterling.o MSorig.o MSorigc.o MStrans.97.o convert.o \
	frausi.o p46merdiff.o g46merdiff.o lj.o modperm.o modf1com.o mododesd.o EAMal.o Alglue.o Mgglue.o \
	Gupta.o orient.o Natb.o sticky.o enumerate.o minperm.o minpermdist.o LB2.o \
	dgetrf.o dgetri.o reorient.o thomson.o Q4.o basinsampling.o tether.o tetherfunc.o BLN.o \
	gauss.o newmindist.o centrecom.o qorderparam_blj.o qorderparam_lj.o bspt.o GMINdump.o quad.o \
	rigidbaa.o checkd.o dbpg.o dbptd.o dmblmorse.o gbcalamitic.o gbdiscotic.o gem.o gbddp.o linrod.o \
        lwotp.o	newcapsid.o newpah.o msgayberne.o mstbin.o multisitepy.o multstock.o paha.o pahw99.o pap.o pyg.o silane.o \
        multpaha.o newtip.o patchy.o asaoos.o pygdp.o stockaa.o tetrahedra.o waterpdc.o waterpkz.o takestepmsgb.o takestepmspy.o \
	gay-berne.o JM.o projI.o projIh.o model1.o FS.o vgw.o \
	mysd.o qdtest.o qdtest2.o MB.o dlsode.o dummyoptim.o Ackland_metals.o Ackland_wrapper.o DF1.o CSMinit.o \
	read_cmd_args.o display_version.o file_manager.o output.o chaperonin.o bulkmindist.o perc.o \
	ga_modules.o ga_main.o ga_select.o ga_bln.o ga_cluster.o

OBJS = ${OBJS1} ${OBJS2}
GENF90FILES = porfuncs.f90
CHDUM = chdummy.o
AMHDUM = amhdummy.o
AMB9DUM = amber9dummy.o
VPATH = .:CHARMM:AMH:.
LDFLAGS = -L.
DEFS =
# CPP = /usr/bin/cpp
CPP = /lib/cpp
CPFLAGS = -traditional -P

# }}}
#
# location of the source for the AMBER9 interface
#
##### AMBER9 ##################### {{{
AMB9SRC = ../../AMBER
SANDERSRC = ../../AMBER/src/sander
# }}}
#
# If you are not building CGMIN.X.X comment everything inside the 31 block!
# You MUST comment the BLAS_EXCLUDE_LIST if compiling GMIN, AMBGMIN and AMHGMIN
#

# Define which CHARMM version you want to use with GMIN. Just uncomment the option you want. 
# for CHARMM 31 (recommended)
#   CTYPE=C31
# for CHARMM 35 (working, but in testing right now)
#   CTYPE=C35
# No CHARMM
   CTYPE=C0
#
# set CHARMM build environment variable and compilation command
# Only used for C31 and C35
# Change PGF90 to ifort to use the ifort compiler etc.
#
chcommand = ./install.com gnu medium keepo keepf PGF90 OPT AMD

ifeq (${CTYPE},C31)
#### start CHARMM 31 {{{
    C31SRC = /home/${USER}/svn/CHARMM31
##  C31SRC = /home/${USER}/svn/charmm31.mpi
    BLAS_EXCLUDE_LIST = dnrm2.o daxpy.o dcopy.o ddot.o
    FCMDIR = ${C31SRC}/source/fcm
    SRC31 = charmm_main.src # energy.src
    EXTRAS = myblas.o mylapack.o
    PREFLX = ${C31SRC}/tool/prefx_gnu
#
# The following line should probably never be uncommented!
##  PREFLX = ${C31SRC}/tool/preflx
    PREFDIR = ${C31SRC}/build/gnu
##### CHARMM 31 block end 
    LIBDIR31=${C31SRC}/lib/gnu
    CHOBJS31 = $(LIBDIR31)/help.o $(LIBDIR31)/iniall.o $(LIBDIR31)/miscom.o $(LIBDIR31)/usersb.o
    CHLIBS31 = $(LIBDIR31)/adumb.a \
	        $(LIBDIR31)/flucq.a $(LIBDIR31)/cadint.a $(LIBDIR31)/cheq.a $(LIBDIR31)/cff.a $(LIBDIR31)/correl.a $(LIBDIR31)/dimb.a \
	        $(LIBDIR31)/emap.a $(LIBDIR31)/dynamc.a $(LIBDIR31)/energy.a $(LIBDIR31)/gamint.a $(LIBDIR31)/gukint.a \
		   $(LIBDIR31)/gener.a $(LIBDIR31)/image.a $(LIBDIR31)/io.a $(LIBDIR31)/machdep.a $(LIBDIR31)/manip.a $(LIBDIR31)/mbond.a \
		   $(LIBDIR31)/mc.a $(LIBDIR31)/minmiz.a $(LIBDIR31)/misc.a $(LIBDIR31)/mmff.a $(LIBDIR31)/molvib.a $(LIBDIR31)/nbonds.a \
		   $(LIBDIR31)/pert.a $(LIBDIR31)/quantum.a $(LIBDIR31)/rxncor.a $(LIBDIR31)/shapes.a $(LIBDIR31)/solvation.a \
		   $(LIBDIR31)/util.a $(LIBDIR31)/vibran.a libcharmm.a
#### end CHARMM 31 }}}
endif

ifeq (${CTYPE},C35)

#### start CHARMM 35 {{{
    C35SRC = /home/${USER}/svn/CHARMM35
    BLAS_EXCLUDE_LIST = dnrm2.o daxpy.o dcopy.o ddot.o
    FCMDIR = ${C35SRC}/source/fcm
    SRC35 = charmm_main.src # energy.src
    EXTRAS = myblas.o mylapack.o
    PREFLX = ${C35SRC}/tool/prefx_gnu
    PREFDIR = ${C35SRC}/build/gnu
    LIBDIR35=${C35SRC}/lib/gnu
    CHOBJS35 = $(LIBDIR35)/help.o $(LIBDIR35)/iniall.o $(LIBDIR35)/miscom.o $(LIBDIR35)/usersb.o
    CHLIBS35 = $(LIBDIR35)/adumb.a \
	        $(LIBDIR35)/flucq.a $(LIBDIR35)/cadint.a $(LIBDIR35)/cheq.a $(LIBDIR35)/cff.a $(LIBDIR35)/correl.a $(LIBDIR35)/dimb.a \
	        $(LIBDIR35)/emap.a $(LIBDIR35)/dynamc.a $(LIBDIR35)/energy.a $(LIBDIR35)/gamint.a $(LIBDIR35)/gukint.a \
		   $(LIBDIR35)/gener.a $(LIBDIR35)/image.a $(LIBDIR35)/io.a $(LIBDIR35)/machdep.a $(LIBDIR35)/manip.a $(LIBDIR35)/mbond.a \
		   $(LIBDIR35)/mc.a $(LIBDIR35)/minmiz.a $(LIBDIR35)/misc.a $(LIBDIR35)/mmff.a $(LIBDIR35)/molvib.a $(LIBDIR35)/nbonds.a \
		   $(LIBDIR35)/pert.a $(LIBDIR35)/quantum.a $(LIBDIR35)/rxncor.a $(LIBDIR35)/shapes.a $(LIBDIR35)/solvation.a \
		   $(LIBDIR35)/util.a $(LIBDIR35)/vibran.a libcharmm.a
#### end CHARMM 35 }}}

endif

###################################### COMPILERS AND COMPILER FLAGS ###################################### {{{
#
######## The Portland Group Compiler Technology Fortran 90 compiler {{
#FC = pgf90
# FC = mpif77  # for lam - don't forget to uncomment MPI!
FC = mpif90  # for mpich and openmpi - don't forget to uncomment MPI!
DEFS = -DMPI
# The usual flags for AMBGMIN:
# FFLAGS= -Mextend -O3 -Mvect=assoc,cachesize:1024000,recog,transform
# These are the CHARMM31 flags for mpif90 64 bit library.
# It is ESSENTIAL to use the same flags as for the CHARMM build!!!!
#
FFLAGS= -Mextend -O3 -Munroll -Mnoframe 
# FFLAGS= -Mextend -O0 -Mnoframe
# FFLAGS= -Mextend -C -g -traceback
# FFLAGS= -Mextend -g -traceback
# Debugging flags
#  FFLAGS= -Mextend -C -g -gopt -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mcoff -Mdwarf1 -Mdwarf2 -Mdwarf3 -Melf -Mnodwarf -Mpgicoff -traceback
# flags for AMBER9
#
FREEFORMAT_FLAG= -Mfree
EXTRA_FLAGS=-module
#
NOOPT = -O0 -Mextend
SEARCH_PATH =  -I..
# the double reference to -lblas seems to be needed!
LIBS = libmyblas.a libmylapack.a libmyblas.a
LDFLAGS= -L.
SWITCH=pgi
# }}}
###### end of The Portland Group Compiler Technology Fortran 90 compiler }}}
#
### NAGWare Fortran {{{ 
# use this compiler for nag/64/5.1
# FC = f95
# DEFS+=-DNAG
# use this compiler for nag/64/5.2
# FC = nagfor
# FFLAGS = -132 -maxcontin=3000 -kind=byte -mismatch_all -O0
# FFLAGS = -132 -maxcontin=3000 -kind=byte -mismatch_all -O3 
# this line is for garden variety debugging 
# FFLAGS = -132 -maxcontin=3000 -C -g -kind=byte -mismatch_all
# this line is for thorough but slow debugging 
# FFLAGS = -132 -maxcontin=3000 -C=all -mtrace=all -gline -kind=byte
# NOOPT= -O0 -132  -kind=byte
# SEARCH_PATH = -I..
# LDFLAGS= -L.
# SWITCH=nag
#
# the double reference to -lblas seems to be needed!
#
# LIBS = libmyblas.a libmylapack.a libmyblas.a 
#
# flags for AMBER9
#
# FREEFORMAT_FLAG= -free
# EXTRA_FLAGS=-I
#
###### end of NAGWare Fortran 95 compiler flags }}}
#
##########################################################
# Intel compilers {{{
#
# GMIN + CHARMM35 - tested with the options which contain "##" in front 
# 
# FC = ifort
# FC = mpif77 
# FC = mpif90  
# DEFS = -DMPI 
##### ifort debugging flags 
# FFLAGS= -132 -C -g -traceback -debug full
# FFLAGS= -132 -O0 -g -traceback -fpe:0 -check all
# FFLAGS= -132 -g -debug all -check all -implicitnone -warn unused -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds
##### ifort optimization flags 
# FFLAGS= -132 -Vaxlib -O3 # for ifc
# FFLAGS= -132 -O3 -ip -static # -ipo or -fast doesnt work 
# FFLAGS= -132 -O4
##FFLAGS= -O2 -extend_source -g -traceback
# NOOPT= -132 -O0
# SWITCH=ifort
# SEARCH_PATH= -I..
# LDFLAGS= -L.
# LIBS = libmyblas.a libmylapack.a libmyblas.a  
# FREEFORMAT_FLAG= -free
# EXTRA_FLAGS=-I
#
# End of intel compilers }}}
### Gfortran  {{{
#
# FC = gfortran
# FC = mpif90
# DEFS = -DMPI
#  FFLAGS= -ffixed-line-length-265 -O0 
#  FFLAGS= -ffixed-line-length-132 -O3 -ftree-vectorize
# FFLAGS= -ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -ftrapv 
# FFLAGS= -ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic
# NOOPT= -O0 -ffixed-line-length-132
# SWITCH=gfortran
# SEARCH_PATH =  -I..
# LIBS = libmyblas.a libmylapack.a libmyblas.a 
# LDFLAGS = -LBLAS -LLAPACK
# FREEFORMAT_FLAG= -ffree-form
# EXTRA_FLAGS=-I
#
#  End Gfortran }}}
# }}}
###############################################################################################
# Pass the subversion revision number into the code to print it in the output
   #DEFS+=-DSVNVERSION="`./svn_revision.sh`"
###################################### RULES AND TARGETS ###################################### {{{
.SUFFIXES:
.SUFFIXES: .o .f .F .f90

.f90.o:
	$(FC) $(FFLAGS) ${SEARCH_PATH} -c $<
.f.o:
	$(FC) $(FFLAGS) ${SEARCH_PATH} -c $<
.F.f:
	$(CPP) $(CPFLAGS) $(DEFS) $< > $@
.F90.f90:
	$(CPP) $(CPFLAGS) $(DEFS) $< > $@

default: $(GPROG)

GMIN: $(GPROG)
gmin: $(GPROG)

CHGMIN: $(CPROG)
chgmin: $(CPROG)


$(GPROG): $(CHDUM) $(AMHDUM) $(AMB9DUM) $(OBJS) $(EXTRAS) 
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ $(EXTRAS) $(OBJS) $(CHDUM) $(AMHDUM) $(AMB9DUM) $(LDFLAGS) $(LIBS)

ifeq (${CTYPE},C31)
$(CPROG): $(AMHDUM) $(AMB9DUM) $(OBJS) $(EXTRAS) libcharmm.a 
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ ${CHOBJS31} $(EXTRAS) $(OBJS) $(AMHDUM) $(AMB9DUM) $(LDFLAGS) $(LIBS) \
        ${CHLIBS31} ${CHLIBS31} ${CHLIBS31}
endif
ifeq (${CTYPE},C35)
$(CPROG): $(AMHDUM) $(AMB9DUM) $(OBJS) $(EXTRAS) libcharmm.a 
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ ${CHOBJS35} $(EXTRAS) $(OBJS) $(AMHDUM) $(AMB9DUM) $(LDFLAGS) $(LIBS) \
        ${CHLIBS35} ${CHLIBS35} ${CHLIBS35}
endif

AMBGMIN: $(AMBPROG) 
ambgmin: $(AMBPROG) 

$(AMBPROG): $(OBJS) $(EXTRAS) $(CHDUM) $(AMHDUM) libamber.a
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ $(OBJS) $(EXTRAS) $(CHDUM) $(AMHDUM) $(LDFLAGS) libamber.a $(LIBS) 

AMHGMIN: $(AMHPROG)
amhgmin: $(AMHPROG)

$(AMHPROG): $(OBJS) $(EXTRAS) $(CHDUM) $(AMB9DUM) libamh.a
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ $(OBJS) $(EXTRAS) $(CHDUM) $(AMB9DUM) $(LDFLAGS) libamh.a $(LIBS) 

#  no agressive optimizations for selected files to cut down on compile time
keyword.o: keyword.f
	${FC} ${NOOPT} ${SEARCH_PATH} -c keyword.f

clean:
	rm -f ${GPROG} ${CPROG} $(AMHPROG) $(AMBPROG) $(OBJS) *.mod $(EXTRAS) *.lst *.o pref.dat prefx.msg nag.f90 FOR021.DAT porfuncs.f90 display_version.f90 mc.f ptmc.f bspt.f main.f *.a
	if test -d ../../BLAS ;  then cd ../../BLAS ; make clean ; fi
	if test -d ../../LAPACK ;  then cd ../../LAPACK ; make clean ; fi
	if test -d CHARMMinterface ;  then cd CHARMMinterface ; make clean ; fi
	if test -d CHARMMinterface35 ;  then cd CHARMMinterface35 ; make clean ; fi
	if test -d AMH ;  then cd AMH ; make clean ; fi
	if test -d $(SANDERSRC) ; then cd $(SANDERSRC) ; make clean ; fi

cleanexe:
	rm -f $(CPROG) $(GPROG) $(AMHPROG) $(AMBPROG)

timing:
	rm -f GMIN.2.0
	$(FC) $(FFLAGS) -p $(OBJS) $(EXTRAS) -o $(GPROG) $(LIBS)

feedback:
	$(FC) $(FFLAGS) -xprofile=use:gmin.profile *.f -o $(GPROG)

rebuild:
	make clean
	make

rebuildamh:
	make clean
	make AMHGMIN

libamber.a:
	export SRCDIR=$(CURDIR);cd ${SANDERSRC}; make lib1 FC="${FC}" FFLAGS="${FFLAGS}" \
	FREEFORMAT_FLAG="${FREEFORMAT_FLAG}" EXTRA_FLAGS="${EXTRA_FLAGS}"
libamh.a: SAT-Ghost
	cd AMH; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}" 
libmyblas.a: SAT-Ghost
	cd ../../BLAS; make double FC="${FC}" FFLAGS="${FFLAGS}" BLAS_EXCLUDE_LIST="${BLAS_EXCLUDE_LIST}";\
	cp libmyblas.a ../GMIN/source/
libmylapack.a: SAT-Ghost
	cd ../../LAPACK; make selection FC="${FC}" FFLAGS="${FFLAGS}" NOOPT="${NOOPT}";\
	cp libmylapack.a ../GMIN/source/


ifeq (${CTYPE},C31)
libcharmm.a: SAT-Ghost
	echo ${chcommand}
	cd ${C31SRC}; ${chcommand}
	echo "Makefile> tail of" ${PREFDIR}/gnu.log
	tail -2 ${PREFDIR}/gnu.log
	cd CHARMMinterface; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}" PREFLX="${PREFLX}" \
	PREFDIR="${PREFDIR}" \
	CTYPE="${CTYPE}" FCMDIR="${FCMDIR}" C31SRC="${C31SRC}" SRC31="${SRC31}"
endif
ifeq (${CTYPE},C35)
libcharmm.a: SAT-Ghost
	echo ${chcommand}
	cd ${C35SRC} ; ${chcommand}
	echo "Makefile> tail of" ${PREFDIR}/gnu.log
	tail -2 ${PREFDIR}/gnu.log
	cd CHARMMinterface35; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}" PREFLX="${PREFLX}" \
	PREFDIR="${PREFDIR}" \
	CTYPE="${CTYPE}" FCMDIR="${FCMDIR}" C35SRC="${C35SRC}" SRC35="${SRC35}"
endif

SAT-Ghost:

porfuncs.f90: porfuncs.csh
	./porfuncs.csh ${SWITCH} > porfuncs.f90

# include ../../include/Makefile.gmin.display_version 

DISPLAY_VERSION=../../SCRIPTS/all/display_version.sh

dvopts=fflags "${FFLAGS}" prog GMIN fc_full_name "${FULL_COMPILER_NAME}" fc_exec "${FC}" make_opts "${MAKE_OPTS}" 

display_version.f90: ${DISPLAY_VERSION} 
	${DISPLAY_VERSION} ${dvopts} > display_version.f90

#read_cmd_args.f90: read_cmd_args.F90

# }}}
###################################### DEPENDENCIES ###################################### {{{
${GPROG}: libmyblas.a libmylapack.a
${CPROG}: libmyblas.a libmylapack.a libcharmm.a
${AMHPROG}: commons.o libmyblas.a libmylapack.a libamh.a
${AMBPROG}: libmyblas.a libmylapack.a libamber.a porfuncs.o
libcharmm.a: commons.o modcharmm.o modmxatms.o
libamber.a: commons.o modamber9.o porfuncs.o 
#${OBJS2}: ${OBJS1} 

Alglue.o:      commons.o
BLJcluster.o:  commons.o
EAMLJ_sub.o:   commons.o
EAMal.o:       commons.o
Gupta.o:       commons.o
MSorig.o:      commons.o modconsts.o dist.o
MSorigc.o:     commons.o modconsts.o dist.o
MStrans.97.o:  commons.o modconsts_trans_97.o dist.o
Mgglue.o:      commons.o
PachecoC60.o:  commons.o
Pbglue.o:      commons.o
SW.o:          commons.o
SiSW.o:        commons.o
Zetterling.o:  commons.o
amber.o:       commons.o modamb.o
bmin.o:        commons.o modamb.o
capsid.o:      commons.o
centre.o:      commons.o
centrecom.o:      commons.o
cgmin.o:       commons.o
checkd.o:      commons.o
compress.o:    commons.o
compress2.o:   commons.o
countatoms.o:  modamber9.o
dbpg.o:        commons.o
dbptd.o:       commons.o
dmblmorse.o:   commons.o
ddfpmin.o:     commons.o
dftb.o:        commons.o
dzugutov.o:    commons.o
evstep.o:      commons.o
f1dim.o:       commons.o modf1com.o
fd.o:          commons.o
fdm.o:         commons.o
finalio.o:     commons.o modamb.o qmod.o modcharmm.o AMH/amhglobals.o modamber9.o
finalq.o:      commons.o qmod.o
frausi.o:      dist.o
gbcalamitic.o: commons.o
gbdiscotic.o:  commons.o
gbddp.o:       commons.o
gem.o:         commons.o
grnd.o:        commons.o
hmat1n.o:      commons.o
hmatd_.o:      commons.o
hsmove.o:      commons.o
io1.o:         commons.o modamb.o modperm.o qmod.o modcharmm.o porfuncs.o
keyword.o:     commons.o modamb.o modcharmm.o porfuncs.o ga_modules.o
linmin.o:      commons.o modf1com.o
linrod.o:      commons.o
ljcoulomb.o:   commons.o
lj.o:          commons.o
ljpbin.o:      commons.o
ljpshift.o:    commons.o
lwotp.o:       commons.o
newpah.o:      commons.o
main.f: main.F
main.o:        commons.o modf1com.o countatoms.o modperm.o qmod.o modamb.o modmxatms.o porfuncs.o 
mcruns.o:      commons.o
mf.o:          commons.o
mnbrak.o:      commons.o
commons.o:     countatoms.o
dummyoptim.o:  commons.o
morse.o:       commons.o
mylbfgs.o:     commons.o modamb.o porfuncs.o
mstbin.o:      commons.o
multstock.o:   commons.o
newtip.o:      commons.o
newcapsid.o:   commons.o
patchy.o:      commons.o
asaoos.o:      commons.o
odesd.o:       commons.o mododesd.o
olami.o:       commons.o
otp.o:         commons.o
BLN.o:         commons.o
pah.o:         commons.o
multpaha.o:    commons.o qmod.o
pahw99.o:      commons.o
paha.o:        commons.o
pap.o:         commons.o
potential.o:   commons.o modperm.o qmod.o modcharmm.o porfuncs.o
ptgrp.o:       commons.o
pyg.o:         commons.o
pygdp.o:       commons.o
quench.o:      commons.o qmod.o porfuncs.o
rad.o:         commons.o
rdpot.o:       commons.o
saveit.o:      commons.o qmod.o
sc.o:          commons.o
seed.o:        commons.o
silane.o:      commons.o
stockaa.o:     commons.o
strand.o:      commons.o
supermc.o:     commons.o
symmetry.o:    commons.o porfuncs.o
symmetrycsm.o:    commons.o porfuncs.o
taboo.o:       commons.o
takestep.o:    commons.o
takestep2.o:   commons.o
tetrahedra.o:  commons.o
tip.o:         commons.o
tosifumi.o:    commons.o
welch.o:       commons.o
waterpdc.o:    commons.o
waterpkz.o:    commons.o
zwischen.o:    commons.o modf1com.o 
stock.o:       commons.o
sticky.o:      commons.o
tether.o:      tetherfunc.o
mycpu_time.o:  commons.o
gauss.o:       commons.o
mc.o:          qmod.o modcharmm.o porfuncs.o mc.f mc.F AMH/amhglobals.o AMH/amh_interfaces.o AMH/E_write.o commons.o operations.o
mc.f: mc.F
ptmc.f: ptmc.F
ptmc.o:          qmod.o modcharmm.o porfuncs.o ptmc.f ptmc.F
bspt.f: bspt.F
bspt.o:          qmod.o modcharmm.o porfuncs.o bspt.f bspt.F
GMINdump.o: commons.o qmod.o porfuncs.o output.o
output.o: commons.o file_manager.o
quad.o:        commons.o
enumerate.o:  commons.o
convert.o:  porfuncs.o
getorbits.o: commons.o
mysd.o: commons.o
qdtest.o: commons.o
MB.o: commons.o 
sort2.o: qmod.o commons.o
model1.o: commons.o
FS.o: commons.o
Ackland_wrapper.o: commons.o
DF1.o: commons.o
CSMinit.o: commons.o
minpermdist.o: commons.o porfuncs.o
grouprotation.o: commons.o
multisitepy.o: gay-berne.o
finalio.o: gay-berne.o
chaperonin.o: commons.o
bulkmindist.o: commons.o
read_cmd_args.o: porfuncs.o display_version.o read_cmd_args.f90 commons.o
ga_main.o: ga_modules.o commons.o
ga_select: ga_modules.o commons.o
ga_bln.o: ga_modules.o commons.o
ga_cluster.o: ga_modules.o commons.o

# }}}
