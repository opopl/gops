#vim:ft=make

PROGNAME=GMIN
GENFFILES = porfuncs.f90 mc.f ptmc.f bspt.f main.f 
NOTUSEDSOURCE= io_sanity.f90 molinfo.f90 transforms.f90
VPATH = .:AMH:.
#
#
PPATH=$(PWD)
INCPATH = $(PPATH)/../INCLUDE/
include $(INCPATH)/make.inc
#
#
include $(INCPATH)/targets
include $(INCPATH)/mkf

$(PROG): $(LIBS)

libamh.a: SAT-Ghost
	cd AMH; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}" 

libamh.a: amhglobals.o

include deps
read_cmd_args.o: display_version.o
