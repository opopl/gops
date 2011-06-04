#vim:ft=make

PROGNAME=GMIN
GENFFILES = porfuncs.f90 dv.f90 
#
#
PPATH=$(PWD)
INCPATH = $(PPATH)/../INCLUDE/
include $(INCPATH)/make.inc
#
#
include $(INCPATH)/targets

$(PROG): $(LIBS)

include $(DEPS)
