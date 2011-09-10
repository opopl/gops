#vim:ft=make

PROGNAME=gmi
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

$(DEPS): rca.f90 dv.f90 porfuncs.f90

include $(DEPS)
