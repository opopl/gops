#vim:ft=make

PROGNAME=GMIN
GENFFILES = porfuncs.f90 display_version.f90 mc.f ptmc.f bspt.f main.f 
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

include deps
