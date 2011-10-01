
default: init dirs $(AUXF) deps $(PROG)

.PHONY: default c clean rmdep deps tg bindir moddir libdir objdir dirs init
.PHONY: pgi nag gf ifort debug opt dx dxdir

#SUFFIXES..., F => o {{{

.SUFFIXES:
.SUFFIXES: .o .f .F .f90

%.o: %.f90
	@echo $< 
	$(FC) $(FFLAGS) $(SEARCH_PATH) -c $< -o $@
	cp $@ $(OBJSPATH)

%.o: %.f
	@echo $< 
	$(FC) $(FFLAGS) $(SEARCH_PATH) -c $< -o $@
	cp $@ $(OBJSPATH)

.F.f:
	$(CPP) $(CPFLAGS) $(DEFS) $< > $@
.F90.f90:
	$(CPP) $(CPFLAGS) $(DEFS) $< > $@

#}}}
#compilers {{{
pgi gf nag ifort: 
	$(MKCOMP) $@ > $(COMP)
debug opt: 
	$(MKCOMP) $@ >> $(COMP)
$(COMP): 
	touch $@

#}}}
# init deps porfuncs rca dv header {{{

init: 
	touch $(INITFILES)

porfuncs.f90: $(PRF)
	$(PRF) $(SWITCH) > $@

#rca.f90: $(RCA)
	#$(RCA) $(PROGNAME) > $@

dv.f90: $(DV)
	$(DV) $(DVOPTS) > $@

header.f90: $(HEADER) 
	$(HEADER) > $@

#}}}
#deps help tg prog bindir cup {{{

deps: $(DEPS)
	$(MKDEP) $@

rmdep:
	rm $(DEPS)

tg:
	ctags -R $(SOURCE)

cup:
	rm -rf $(NOTUSEDSOURCE)

bindir: 
	mkdir -p $(BINPATH)
moddir: 
	mkdir -p $(MODPATH)
libdir: 
	mkdir -p $(LIBAPATH)
objdir: 
	mkdir -p $(OBJSPATH)

dirs: bindir moddir libdir objdir
#}}}
# PROGNAME* {{{

# o => optimization
# g => debug
opts:= o g

$(PROGNAME)_o: dirs $(DEPS) $(OBJS) 
	$(FC) $(FFLAGS_o)  -o $@ $(OBJS) $(LDFLAGS) $(LIBS)
	cp $(PROG) ./
	cp $(PROG) $(PROG)_o_$(FC)

$(PROGNAME)_g: dirs $(DEPS) $(OBJS) 
	$(FC) $(FFLAGS_g)  -o $@ $(OBJS) $(LDFLAGS) $(LIBS)
	cp $(PROG) ./
	cp $(PROG) $(PROG)_g_$(FC)

$(PROGNAME): $(PROG)
$(PROGNAME)_use_base: dirs $(OBJS) 
	$(FC) $(FFLAGS)  -o $@ $(LPU_OBJS) $(LDFLAGS) $(LPBASE) $(LIBS)
	cp $@ $(BINPATH)

# }}}
# DEPS PROG {{{

$(DEPS): $(MKDEP)
	$(MKDEP) $@

$(PROG): dirs $(DEPS) $(OBJS) 
	rm -f $(PROG) dv.f90
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)
	cp $(PROG) ./
	cp $(PROG) $(PROG)_$(FC)

#}}}
#libs {{{

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
#}}}
dx: $(dxfile) dxdir
	$(DXSH)	$(PROGNAME) > $(dxfile)

#clean cleanall{{{

cbase:
	rm -f $(LPBASE)

c: clean

clean:
	rm -f $(PROG)
	rm -f $(DEPS) $(GENFFILES)
	rm -f $(OBJS) $(AUXF) $(GENFFILES) *.mod $(MODPATH)/*.mod *.mod *.oo *.o *.ipa
 
cleanall:
	rm -f $(PROG) $(OBJS) ${GENFFILES} $(LLIBS) *.o *.a *.lst
	if test -d $(BLASPATH) ;  then cd $(BLASPATH) ; make clean ; fi
	if test -d $(LAPACKPATH) ;  then cd $(LAPACKPATH) ; make clean ; fi
	if test -d NEB ;  then cd NEB ; make clean ; fi
	if test -d CONNECT ;  then cd CONNECT ; make clean ; fi
	if test -d AMH ;  then cd AMH ; make clean ; fi

#}}}
