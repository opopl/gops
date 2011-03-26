#vim:ft=make
# Definitions {{{
#

ROOT_PATH=$(PWD)/../
INC_PATH=$(ROOT_PATH)/INCLUDE/

#PROGNAME=$(shell sh -c "readlink -f . | awk -F '/' '{ print $NF }'" )
PROGNAME=GMIN
GENFFILES = porfuncs.f90 display_version.f90 mc.f ptmc.f bspt.f main.f 

include $(INC_PATH)/def

AMHPROG = $(BIN_PATH)/AMHGMIN

# }}}
# Object files {{{

ALLSOURCE := $(wildcard *.f) $(wildcard *.f90) $(wildcard *.F)

## NOTUSEDSOURCE: source files not used in present GMIN build {{{
#NOTUSEDSOURCE=\
#blas.f\
#bmin.f\
#bsptsave.F\
#charmmBildc.f\
#compress2.f\
#dblpy.f90\
#dfuncs.f90\
#dgetf2.f\
#diffract.f90\
#dlaswp.f\
#dtrti2.f\
#dtrtri.f\
#evstep.f \
#fsc1TH.f90\
#funcs.f90\
#io.f90\
#lapack.f\
#lbfgs.f\
#mek-quake.f90\
#mf.f\
#mind.f\
#nonATLASblas.f\
#olami.f \
#otp.f\
#ptbs.f\
#rdpot.f \
#rk45ad.f\
#symmetry3.f90\
#takestep2.f\
#tetherfuncall.f90 \
#wenzel.1D.f\
#wenzel.2D.f\
#wham_djw.f90\
## }}}

# NOTUSEDSOURCE: source files not used in present GMIN build {{{
NOTUSEDSOURCE= io_sanity.f90 molinfo.f90 transforms.f90
## }}} 

SOURCE=$(filter-out $(NOTUSEDSOURCE),$(ALLSOURCE))
OBJS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCE))))

# object files from source files generated during compilation
OBJSrt := porfuncs.o display_version.o
OBJS+= $(OBJSrt)

AMHDUM = amhdummy.o
VPATH = .:AMH:.

# }}}

include $(INC_PATH)/targets
include $(INC_PATH)/mkf

## DEPENDENCIES  {{{
$(PROG): $(LIBS)
${AMHPROG}: commons.o $(LIBS) libamh.a

Alglue.o:      commons.o
BLJcluster.o:  commons.o
EAMLJ_sub.o:   commons.o
EAMal.o:       commons.o
Gupta.o:       commons.o
MSorig.o:      commons.o modconsts.o dist.o
MSorigc.o:     commons.o modconsts.o dist.o
MStrans.97.o:  commons.o modconsts_trans_97.o dist.o
Mgglue.o:      commons.o
PachecoC60.o:  commons.o
ptbs.o: commons.o tetherfunc.o 
Pbglue.o:      commons.o
SW.o:          commons.o
SiSW.o:        commons.o
Zetterling.o:  commons.o
amber.o:       commons.o modamb.o
capsid.o:      commons.o
centre.o:      commons.o
centrecom.o:      commons.o
cgmin.o:       commons.o
checkd.o:      commons.o
compress.o:    commons.o
countatoms.o:  modamber9.o
display_version.o: ../../SCRIPTS/all/display_version.sh
dbpg.o:        commons.o
dbptd.o:       commons.o
dmblmorse.o:   commons.o
ddfpmin.o:     commons.o
dftb.o:        commons.o
dzugutov.o:    commons.o
evstep.o:      commons.o
f1dim.o:       commons.o modf1com.o
fd.o:          commons.o
fdm.o:         commons.o
finalio.o:     commons.o modamb.o qmod.o modcharmm.o AMH/amhglobals.o pymodule.o
finalq.o:      commons.o qmod.o
frausi.o:      dist.o
gbcalamitic.o: commons.o
gbdiscotic.o:  commons.o
gbddp.o:       commons.o
gem.o:         commons.o
grnd.o:        commons.o
hmat1n.o:      commons.o
hmatd_.o:      commons.o
hsmove.o:      commons.o
io1.o:         commons.o modamb.o modperm.o qmod.o modcharmm.o porfuncs.o
keyword.o:     commons.o modamb.o modcharmm.o porfuncs.o
linmin.o:      commons.o modf1com.o
linrod.o:      commons.o
ljcoulomb.o:   commons.o
lj.o:          commons.o
ljpbin.o:      commons.o
ljpshift.o:    commons.o
lwotp.o:       commons.o
newpah.o:      commons.o
main.o:        commons.o modf1com.o countatoms.o modperm.o qmod.o modamb.o modmxatms.o porfuncs.o 
mcruns.o:      commons.o
mf.o:          commons.o
mnbrak.o:      commons.o
commons.o:     countatoms.o modmxatms.o
dummyoptim.o:  commons.o
morse.o:       commons.o
mylbfgs.o:     commons.o modamb.o porfuncs.o
mstbin.o:      commons.o
multstock.o:   commons.o
newtip.o:      commons.o
odesd.o:       commons.o mododesd.o
olami.o:       commons.o
otp.o:         commons.o
BLN.o:         commons.o
pah.o:         commons.o
multpaha.o:    commons.o qmod.o
pahw99.o:      commons.o
paha.o:        commons.o
potential.o:   commons.o modperm.o qmod.o modcharmm.o porfuncs.o
ptgrp.o:       commons.o
pyg.o:         commons.o
pygdp.o:       commons.o
quench.o:      commons.o qmod.o porfuncs.o
rad.o:         commons.o
rdpot.o:       commons.o
saveit.o:      commons.o qmod.o
sc.o:          commons.o
seed.o:        commons.o
stockaa.o:     commons.o
strand.o:      commons.o
supermc.o:     commons.o
symmetry.o:    commons.o porfuncs.o
symmetrycsm.o:    commons.o porfuncs.o
taboo.o:       commons.o
takestep.o:    commons.o
takestep2.o:   commons.o
tetrahedra.o:  commons.o
g46merdiff.o: modhess.o
tip.o:         commons.o
tosifumi.o:    commons.o
welch.o:       commons.o
waterpdc.o:    commons.o
waterpkz.o:    commons.o
zwischen.o:    commons.o modf1com.o 
stock.o:       commons.o
sticky.o:      commons.o
tether.o:      tetherfunc.o
mycpu_time.o:  commons.o
gauss.o:       commons.o
mc.o:          qmod.o modcharmm.o porfuncs.o mc.f mc.F AMH/amhglobals.o AMH/amh_interfaces.o AMH/E_write.o commons.o operations.o
bspt.f: bspt.F
main.f: main.F
mc.f: mc.F
ptmc.f: ptmc.F
ptmc.o:          qmod.o modcharmm.o porfuncs.o
bspt.o:          tetherfunc.o commons.o modcharmm.o
GMINdump.o: commons.o qmod.o porfuncs.o
quad.o:        commons.o
enumerate.o:  commons.o
convert.o:  porfuncs.o
getorbits.o: commons.o
mysd.o: commons.o
qdtest.o: commons.o
MB.o: commons.o 
sort2.o: qmod.o commons.o
model1.o: commons.o
FS.o: commons.o
Ackland_wrapper.o: commons.o
DF1.o: commons.o
CSMinit.o: commons.o
minpermdist.o: commons.o porfuncs.o
grouprotation.o: commons.o

# op226 
read_cmd_args.o: porfuncs.o display_version.o commons.o

# }}}
