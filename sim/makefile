#vim:ft=make

PROGNAME=gmin
GENFFILES = porfuncs.f90 dv.f90 
#
#
PPATH=$(PWD)
INCPATH = $(PPATH)/../include/
include $(INCPATH)/make.inc
#
#
include $(INCPATH)/targets

$(PROG): $(LIBS)

include $(DEPS)
