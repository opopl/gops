#vim:ft=make
# Definitions {{{
#

PROGNAME=GMIN
GENFFILES = porfuncs.f90 display_version.f90 mc.f ptmc.f bspt.f main.f 
NOTUSEDSOURCE= io_sanity.f90 molinfo.f90 transforms.f90

VPATH = .:AMH:.

include $(INC_PATH)/targets
include $(INC_PATH)/mkf

$(PROG): $(LIBS)


