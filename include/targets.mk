
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

default: $(AUXF) $(DEPS) $(PROG)

porfuncs.f90: $(PRF)
	$(PRF) $(SWITCH) > $@

rca.f90: $(RCA)
	$(RCA) $(PROGNAME) > $@

dv.f90: $(DV)
	$(DV) $(DVOPTS) > $@

header.f90: $(HEADER) 
	$(HEADER) > $@

$(DEPS): $(MKDEP)
	$(MKDEP) $@

tg:
	ctags -R *.f *.f90 *.F

$(PROGNAME): $(PROG)
	cp $(PROG) ./

bindir: 
	mkdir -p $(BINPATH)

$(PROG): $(OBJS) bindir
	$(FC) $(FFLAGS) $(SEARCH_PATH) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

cup:
	rm -rf $(NOTUSEDSOURCE)

libmyblas.a: SAT-Ghost
	cd $(BLASPATH); make double FC="${FC}" FFLAGS="${FFLAGS}" BLAS_EXCLUDE_LIST="${BLAS_EXCLUDE_LIST}";\
	cp $@ $(PPATH)

libmylapack.a: SAT-Ghost
	cd $(LAPACKPATH); make selection FC="${FC}" FFLAGS="${FFLAGS}" NOOPT="${NOOPT}";\
	cp $@ $(PPATH)

SAT-Ghost:

clean:
	rm -f $(OBJS) $(AUXF) $(GENFFILES) *.mod $(DEPS)
 
cleanall:
	rm -f $(PROG) $(OBJS) ${GENFFILES} *.o *.a *.mod *.lst
	if test -d $(BLASPATH) ;  then cd $(BLASPATH) ; make clean ; fi
	if test -d $(LAPACKPATH) ;  then cd $(LAPACKPATH) ; make clean ; fi
	if test -d NEB ;  then cd NEB ; make clean ; fi
	if test -d CONNECT ;  then cd CONNECT ; make clean ; fi
	if test -d AMH ;  then cd AMH ; make clean ; fi
