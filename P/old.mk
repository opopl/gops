PROG =	../bin/PATHSAMPLE

OBJS =	common.o nodes.o porfuncs.o utils.o key.o dock.o main.o keywords.o setup.o input.o \
	KMCcommit.o mysystem.o sort.o KMC.o inertia.o eig.o amhdump.o charmmdump.o \
	cycle.o mywait.o submitoptimjob.o connectodata.o mindist.o unresdump.o \
	centre.o getnewpath.o addperm.o tssearch.o \
	mindouble.o tsdouble.o Pfold.o NGT.o GT.o getallpaths.o cycle2.o dprand.o sdprnd.o \
        Dijkstra.o getdpair.o Dijinit.o calcorder.o getrpair.o connectd.o getspair.o \
        getpair.o pairdouble.o donedouble.o newmindist.o mergedb.o getupair.o \
	GT2.o GT2data.o GT2FibonacciHeap.o GT2DLL.o GT2input.o GT2FreeMemory.o \
        minperm.o minpermdist.o rigidbodymod.o mathsconstants.o quaternionmatch.o \
        regroupfree.o regroup.o dsort.o savestate.o orderodata.o regroupfree2.o \
	probacc.o newconn.o getfreepair.o getfreebarrier.o kshortestpaths.o Dijinitfly.o \
	getallmin.o myorient.o getusepair.o NGTmem.o NGTrealloc.o NGTrenorm.o \
	NGTremoveid.o NGTremovei.o swapnode.o mymerge.o rbinertia.o rigidb.o diagonalise2.o reweight.o \
	rbperm.o minpermdistrbcom.o Cv.o DOS.o bulkmindist.o frictionfac.o make_conpot.o setup_sis.o rateconst_setup.o \
	checkTS.o remove_unconnected.o  \
	NGT_crstorage.o NGTremovei_crstorage.o NGTremoveid_crstorage.o NGTrenorm_crstorage.o NGTrealloc_crstorage.o

# note that nag64/5.1-216 fails for large memory, but version 365 works
# large memory for pgi and ifort requires the -mcmodel flag!

# WARNING - points.min and point.ts created with ifort executables cannot be
# read by NAG or PGI executables, and vice versa, unless -assume byterecl is used

# Preprocessing
   DEFS =
   CPP = /lib/cpp
   CPFLAGS = -traditional -P
  
# NAG compiler {{{
#  use this compiler for nag/64/5.1
#  FC = f95
#  use this compiler for nag/64/5.2
#FC = nagfor
#  this line is for garden variety debugging 
#  FFLAGS = -132 -maxcontin=3000 -C -g -kind=byte -mismatch_all -ieee=full
#  this line is for thorough but slow debugging 
#  FFLAGS = -132 -maxcontin=3000 -C=all -mtrace=all -gline -kind=byte
# FFLAGS = -132 -maxcontin=3000 -mismatch_all -kind=byte -O0
#FFLAGS = -132 -maxcontin=3000 -mismatch_all -kind=byte -O3 -ieee=full
#NOOPT= -132 -maxcontin=3000 -ieee=full  -kind=byte -O0
#NAGSWITCH = nag
#LIBS = libmyblas.a libmylapack.a libmyblas.a 
#LDFLAGS = -LBLAS -LLAPACK
# }}}
# pathscale compiler {{{
#  FC = pathf95
#  FFLAGS = -extend-source  -g -C
#  FFLAGS =  -extend-source -O3 
#  NOOPT= -extend-source -O0
#  NAGSWITCH = pathscale
#  LIBS = libmyblas.a libmylapack.a libmyblas.a 
#  LDFLAGS = -LBLAS -LLAPACK
# }}}
### Intel  {{{
### please use NAG or PGI if possible. Unformatted points files generated with ifort are incompatible with NAG and PGI
### without the -assume byterecl flag for ifort.
#
#  Without -heap-arrays ifort executables now generate a SEGV in Dijkstra.f90
#
#  FC = ifort
#  FFLAGS= -132 -C -g -heap-arrays -assume byterecl
#  FFLAGS= -132 -g -debug all -check all -implicitnone -warn unused -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds  -assume byterecl
#  FFLAGS= -132 -heap-arrays -O4 -assume byterecl
#  FFLAGS= -132 -O0 -heap-arrays -assume byterecl
#  NOOPT= -O0  -assume byterecl # -mcmodel=large
#  NAGSWITCH=ifort
#  SEARCH_PATH =  -I..
#  LIBS = libmyblas.a libmylapack.a libmyblas.a 
#  LDFLAGS = -LBLAS -LLAPACK
# }}}
##### The Portland Group Compiler Technology Fortran 90 compiler {{{
#
  FC = pgf90
 FFLAGS= -Mextend -O3 -Munroll -Mnoframe
# FFLAGS=  -Mextend -O3 -Munroll -Mnoframe -Mcache_align -Mflushz -mcmodel=medium # -C -g
# FFLAGS= -Mextend -C -g -mcmodel=medium
# NOOPT = -O0 -Mextend -mcmodel=medium
# NAGSWITCH=pgi
 LIBS = -lblas -llapack -lblas
 LDFLAGS = -LBLAS -LLAPACK
#
##### end of The Portland Group Compiler Technology Fortran 90 compiler }}}
### Gfortran  {{{
#
# FC = gfortran
# FFLAGS= -ffixed-line-length-132 -O0 
# FFLAGS= -ffixed-line-length-132 -O2
# FFLAGS= -ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -ftrapv 
# FFLAGS= -ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic
# NOOPT= -O0 -ffixed-line-length-132
# NAGSWITCH=gfortran
# SEARCH_PATH =  -I..
# LIBS = libmyblas.a libmylapack.a libmyblas.a 
# LDFLAGS = -LBLAS -LLAPACK

#
# }}}
###############################################################################################
# Pass the subversion revision number into the code to print it in the output
   DEFS+=-DSVNVERSION="`./version.sh`"
####################################### RULES AND TARGETS ###################################### {{{

# END OF COMPILER SPECIFIC STUFF

.SUFFIXES:
.SUFFIXES: .o .f .F .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<
.c.o:
	$(CC) -c $<
.F.f:
	$(CPP) $(CPFLAGS) $(DEFS) $< > $@
		
# first target encountered is the default build

default: $(PROG) 

$(PROG): $(OBJS) blas_lapack
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

blas_lapack: libmyblas.a libmylapack.a
libmyblas.a:
	cd ../../BLAS; make double FC="${FC}" FFLAGS="${FFLAGS}" BLAS_EXCLUDE_LIST="${BLAS_EXCLUDE_LIST}";\
        cp libmyblas.a ../PATHSAMPLE/source

libmylapack.a:
	cd ../../LAPACK; make selection FC="${FC}" FFLAGS="${FFLAGS}" NOOPT="${NOOPT}";\
        cp libmylapack.a ../PATHSAMPLE/source

porfuncs.f90: porfuncs.csh
	./porfuncs.csh ${NAGSWITCH} > porfuncs.f90

clean:
	rm -f $(PROG) $(OBJS) *.mod porfuncs.f90 main.f *.lst *.a
	cd ../../BLAS; make clean
	cd ../../LAPACK; make clean
$(OBJS): libmyblas.a libmylapack.a
KMC.o:     porfuncs.o common.o
KMCcommit.o:      common.o porfuncs.o
addperm.o:      porfuncs.o common.o
centre.o:      common.o
charmmdump.o:      key.o common.o utils.o
connectodata.o:    key.o common.o
dock.o:         porfuncs.o key.o common.o nodes.o
getallpaths.o:     porfuncs.o common.o key.o utils.o
getallmin.o:     porfuncs.o common.o key.o
getnewpath.o:      porfuncs.o common.o key.o
inertia.o:         common.o rigidbodymod.o
keywords.o:      porfuncs.o nodes.o key.o common.o rigidbodymod.o
main.f: main.F 
main.o:      common.o porfuncs.o rigidbodymod.o dock.o
nodes.o:     porfuncs.o common.o
setup.o:      porfuncs.o utils.o key.o common.o
setup_sis.o:	porfuncs.o utils.o key.o common.o
tssearch.o:   porfuncs.o key.o common.o
unresdump.o:  common.o
Pfold.o:    common.o porfuncs.o
cycle.o:    common.o porfuncs.o
cycle2.o:    common.o porfuncs.o
mind.o: porfuncs.o
mindouble.o: common.o
mysystem.o:   porfuncs.o
mywait.o:     porfuncs.o
nodes.o:           porfuncs.o
submitoptimjob.o:  porfuncs.o nodes.o key.o common.o
tsdouble.o:   common.o
getdpair.o:   common.o
diagonalise2.o: porfuncs.o common.o
Dijkstra.o: common.o porfuncs.o
Dijinit.o: common.o porfuncs.o
Dijinitfly.o: common.o porfuncs.o
connectd.o: common.o
frictionfac.o: common.o
getrpair.o: common.o
getspair.o: common.o porfuncs.o
getupair.o: common.o porfuncs.o
getfreepair.o: common.o porfuncs.o
getpair.o: common.o porfuncs.o
pairdouble.o: common.o 
donedouble.o: common.o 
mergedb.o:     porfuncs.o common.o key.o
GT.o: porfuncs.o common.o GT2.o GT2input.o GT2FreeMemory.o savestate.o
NGT.o: porfuncs.o common.o savestate.o NGTmem.o
NGTrealloc.o: NGTmem.o
NGTrenorm.o: NGTmem.o
NGTremovei.o: common.o NGTmem.o porfuncs.o
NGTremoveid.o: common.o NGTmem.o
NGT_crstorage.o: porfuncs.o common.o savestate.o NGTmem.o
NGTrealloc_crstorage.o: common.o NGTmem.o
NGTrenorm_crstorage.o: NGTmem.o
NGTremovei_crstorage.o: common.o NGTmem.o porfuncs.o
NGTremoveid_crstorage.o: common.o NGTmem.o
GT2.o: common.o GT2data.o GT2FibonacciHeap.o GT2DLL.o
GT2FibonacciHeap.o: GT2data.o
GT2DLL.o: GT2data.o
GT2input.o: GT2data.o GT2DLL.o
GT2FreeMemory.o: GT2data.o GT2DLL.o
minpermdist.o: common.o
minpermdistrbcom.o: common.o
newmindist.o:   common.o rigidbodymod.o
rbinertia.o: common.o
rbperm.o: common.o key.o
rigidb.o: common.o
rigidbodymod.o: mathsconstants.o
regroupfree.o: common.o savestate.o
regroupfree2.o: common.o savestate.o
regroup.o: common.o
orderodata.o: common.o
kshortestpaths: common.o porfuncs.o
reweight.o: porfuncs.o common.o
Cv.o: common.o
DOS.o: common.o
bulkmindist.o: common.o
make_conpot.o: common.o
rateconst_setup.o: common.o porfuncs.o
checkTS.o: common.o porfuncs.o
# END RULES AND TARGETS }}}
