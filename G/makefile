#vim:ft=make

PROGNAME=GMIN
GENFFILES = porfuncs.f90 dv.f90 
#
#
PPATH=$(PWD)
INCPATH = $(PPATH)/../include/

include $(INCPATH)/inc.mk
#
#
include $(INCPATH)/targets.mk

$(PROG): $(LIBS)

$(DEPS): rca.f90 dv.f90 porfuncs.f90 nu.mk

include $(DEPS)
