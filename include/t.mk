
default: dirs $(AUXF) $(DEPS) $(PROG)

#SUFFIXES..., F => o {{{

.SUFFIXES:
.SUFFIXES: .o .f .F .f90

.f90.o:
	$(FC) $(FFLAGS) $(SEARCH_PATH) -c $<
.f.o:
	$(FC) $(FFLAGS) $(SEARCH_PATH) -c $<
.F.f:
	$(CPP) $(CPFLAGS) $(DEFS) $< > $@
.F90.f90:
	$(CPP) $(CPFLAGS) $(DEFS) $< > $@

#}}}
# porfuncs rca dv header {{{

porfuncs.f90: $(PRF)
	$(PRF) $(SWITCH) > $@

rca.f90: $(RCA)
	$(RCA) $(PROGNAME) > $@

dv.f90: $(DV)
	$(DV) $(DVOPTS) > $@

header.f90: $(HEADER) 
	$(HEADER) > $@

#}}}

#deps tg prog bindir cup {{{

deps: rmdep $(DEPS)
rmdep:
	rm $(DEPS)

$(DEPS): $(MKDEP)
	$(MKDEP) $@

tg:
	ctags -R *.f *.f90 *.F

$(PROGNAME): $(PROG)

bindir: 
	mkdir -p $(BINPATH)
moddir: 
	mkdir -p $(MODPATH)
libdir: 
	mkdir -p $(LIBAPATH)

dirs: bindir moddir libdir

b$(PROGNAME): dirs $(OBJS) 
	$(FC) $(FFLAGS) $(SEARCH_PATH) -o $@ $(LPU_OBJS) $(LDFLAGS) $(LPBASE) $(LIBS)
	cp $@ $(BINPATH)

$(PROG): dirs $(OBJS)  
	$(FC) $(FFLAGS) $(SEARCH_PATH) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)
	cp $(PROG) ./

cup:
	rm -rf $(NOTUSEDSOURCE)
#}}}

libamh.a: $(AMH_OBJS)
	cd AMH; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}" 

#blas lapack {{{

blas: $(LBLAS)
lapack: $(LLAPACK)
base: $(LPBASE)

$(LBLAS): $(BLAS_OBJS)
	cd $(BLASPATH); make double FC="${FC}" FFLAGS="${FFLAGS}" BLAS_EXCLUDE_LIST="${BLAS_EXCLUDE_LIST}";\
	mkdir -p $(LIBAPATH)
	cp $(BLASPATH)/libmyblas.a $@

$(LLAPACK): $(LAPACK_OBJS)
	cd $(LAPACKPATH); make selection FC="${FC}" FFLAGS="${FFLAGS}" NOOPT="${NOOPT}";\
	mkdir -p $(LIBAPATH)
	cp $(LAPACKPATH)/libmylapack.a $@

$(LPBASE): $(LPBASE_OBJS) 
	mkdir -p $(LIBAPATH)
	$(ARCH) $(ARCHFLAGS) $(LPBASE) $(LPBASE_OBJS) $^
	$(RANLIB) $@ 

SAT-Ghost:
#}}}

#clean cleanall{{{

cbase:
	rm -f $(LPBASE)

clean:
	rm -f $(OBJS) $(AUXF) $(GENFFILES) $(MODPATH)/*.mod $(DEPS)
 
cleanall:
	rm -f $(PROG) $(OBJS) ${GENFFILES} $(LLIBS) *.o *.a *.lst
	if test -d $(BLASPATH) ;  then cd $(BLASPATH) ; make clean ; fi
	if test -d $(LAPACKPATH) ;  then cd $(LAPACKPATH) ; make clean ; fi
	if test -d NEB ;  then cd NEB ; make clean ; fi
	if test -d CONNECT ;  then cd CONNECT ; make clean ; fi
	if test -d AMH ;  then cd AMH ; make clean ; fi

#}}}
