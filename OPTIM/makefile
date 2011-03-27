
# GPL License Info {{{
# 
#   OPTIM: A program for optimizing geometries and calculating reaction pathways
#   Copyright (C) 1999-2011 David J. Wales
#   This file is part of OPTIM.
#
#   OPTIM is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   OPTIM is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# }}}

# Definitions {{{
#

ROOT_PATH=$(PWD)/../
INC_PATH=$(ROOT_PATH)/INCLUDE/

PROGNAME=OPTIM
GENFFILES = porfuncs.f90 header.f90

include $(INC_PATH)/def

# }}}

# Object files {{{

ALLSOURCE := $(wildcard *.f) $(wildcard *.f90) $(wildcard *.F)
NOTUSEDSOURCE=

SOURCE=$(filter-out $(NOTUSEDSOURCE),$(ALLSOURCE))
OBJS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCE))))

# object files from source files generated during compilation

OBJSrt := porfuncs.o display_version.o
OBJS+= $(OBJSrt)

AMHDUM = amhdummy.o
VPATH = .:AMH:NEB:CONNECT:.

# }}}

LDFLAGS = -L.

#LIBS = -lnn  -lnc libmyblas.a libmylapack.a libmyblas.a libmyblas.a -lnc
LIBS = -lnn  -lnc $(LLIBS) -lnc

include $(INC_PATH)/targets
include $(INC_PATH)/mkf

libnn.a: SAT-Ghost
	cd NEB; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}"
libnc.a: SAT-Ghost
	cd CONNECT; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}"
libamh.a: SAT-Ghost
	cd AMH; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}" 

$(PROG): libnn.a libnc.a $(LLIBS)

libnn.a: modcharmm.o key.o commons.o modunres.o porfuncs.o modmec.o modguess.o efol.o growstring.o gsdata.o intcommons.o intcoords.o amhglobals.o 
libnc.a: libnn.a key.o syminf.o modhess.o amhglobals.o
libamh.a: amhglobals.o

OPTIM.o: commons.o intcommons.o intcoords.o key.o modcharmm.o modguess.o \
	modhess.o modmec.o modneb.o modtwoend.o modunres.o porfuncs.o vecck.o \
	zwk.o
2Dfunc.o: commons.o modhess.o
Ackland_wrapper.o: key.o
adm.o: commons.o key.o porfuncs.o
amber.o: commons.o key.o modamber.o modamber2.o modhess.o
axdiff.o: modhess.o
aziz.o: modhess.o
beig.o: commons.o key.o modneb.o modtwoend.o porfuncs.o
bfgsts.o: commons.o key.o modcharmm.o modhess.o modtwoend.o porfuncs.o \
	vecck.o zwk.o
BLN.o: key.o modhess.o
c10.o: modhess.o
c60diff.o: modhess.o
c60p.o: modhess.o
caardiff.o: modhess.o
capsid.o: modhess.o
changep.o: commons.o key.o modamber.o modneb.o modtwoend.o porfuncs.o
charmmBildc.o: commons.o modamber9.o
Clatmin.o: commons.o
connect.o: commons.o efol.o key.o modamber9.o modcharmm.o modhess.o modneb.o \
	modtwoend.o modunres.o porfuncs.o syminf.o
cpmdlatmin.o: key.o porfuncs.o
ctest.o: modhess.o
dctrap.o: key.o modhess.o
DF1.o: key.o modhess.o
dftb.o: key.o modhess.o
diis.o: gdiis.o key.o
dimer.o: commons.o modhess.o
double.o: modhess.o
dtrap.o: modhess.o
dumpit.o: commons.o key.o
dumpp.o: commons.o key.o modcharmm.o porfuncs.o
dzugutov.o: modhess.o
ectrap.o: key.o
escp.o: modhess.o
EYtrap.o: key.o modhess.o
fd.o: key.o modhess.o
fetchz.o: SD.o amhglobals.o binaryio.o commons.o key.o modamber9.o \
	modcharmm.o modguess.o modmec.o modneb.o modtwoend.o modunres.o \
	msevb_common.o porfuncs.o syminf.o
g2special.o: key.o modhess.o porfuncs.o
g46merdiff.o: modhess.o
geopt.o: amhglobals.o binaryio.o commons.o efol.o key.o modcharmm.o modhess.o \
	modneb.o modtwoend.o porfuncs.o
gmetry.o: commons.o key.o
GoOptim.o: key.o
greatcirc.o: commons.o key.o modcharmm.o modhess.o modtwoend.o modunres.o \
	porfuncs.o zwk.o
h2o.o: modhess.o
hessout.o: modhess.o
hybridmin.o: commons.o key.o modcharmm.o modhess.o porfuncs.o vecck.o zwk.o
inertia.o: commons.o key.o
input.o: porfuncs.o
intbeig.o: commons.o key.o modneb.o modtwoend.o
intbfgsts.o: commons.o key.o modhess.o modunres.o vecck.o zwk.o
intsecdiag.o: commons.o key.o
intxmylbfgs.o: commons.o key.o modtwoend.o
ions.o: modhess.o
JM.o: modhess.o porfuncs.o
keyword.o: amhglobals.o binaryio.o commons.o cubsplstring.o gsdata.o \
	intcommons.o key.o modamber.o modamber9.o modcharmm.o modguess.o \
	modmec.o modmxatms.o modneb.o modtwoend.o modunres.o msevb_common.o \
	porfuncs.o pymodule.o wc.o
ljdiff.o: modhess.o
ljms.o: modhess.o
ljpbin.o: key.o modhess.o
ljpdiff.o: key.o modhess.o
ljpkob.o: key.o modhess.o
ljpshift.o: key.o modhess.o
main.o: commons.o key.o modamber9.o porfuncs.o
MB.o: commons.o modhess.o
mdiff.o: modhess.o
mied.o: modhess.o
mindist.o: key.o porfuncs.o
morph.o: commons.o key.o modcharmm.o modtwoend.o porfuncs.o zwk.o
morse.o: modhess.o
mpdiff.o: key.o modhess.o
msdiff.o: modhess.o
mylbfgs.o: commons.o key.o modamber9.o modcharmm.o modhess.o modtwoend.o \
	modunres.o porfuncs.o specialfuncts.o zwk.o
myorient.o: key.o
odesd.o: commons.o key.o modhess.o porfuncs.o
oldneb.o: commons.o key.o modcharmm.o modhess.o modneb.o modtwoend.o \
	porfuncs.o
p46merdiff.o: modhess.o
p4diff.o: modhess.o
PachecoC60.o: modhess.o
path.o: amhglobals.o commons.o efol.o key.o modcharmm.o modhess.o modunres.o \
	porfuncs.o syminf.o
pertable.o: commons.o key.o
potential.o: SD.o VRTMCY5f.o commons.o finite_differences.o key.o modamber9.o \
	modcharmm.o modhess.o porfuncs.o
projh.o: modhess.o
PV.o: commons.o key.o
rad.o: commons.o key.o
repelsp.o: commons.o key.o
rotd.o: modhess.o
scdiff.o: modhess.o
secdiag.o: commons.o key.o modcharmm.o porfuncs.o
shifth.o: key.o modhess.o
SiO2.o: modhess.o porfuncs.o
sortxyz.o: commons.o
stock.o: commons.o key.o modhess.o vecck.o
SW.o: key.o modhess.o
symmetry.o: commons.o key.o syminf.o
tcheck.o: commons.o
thomson.o: commons.o key.o modhess.o
tight.o: key.o modhess.o
tip.o: commons.o key.o modhess.o
TIPmodes.o: modhess.o
tosifumi.o: modhess.o porfuncs.o
twoend.o: commons.o key.o modcharmm.o modhess.o modtwoend.o porfuncs.o zwk.o
unrescalcdihe.o: commons.o modunres.o
unresconnectsections.o: commons.o key.o modtwoend.o modunres.o
unresdummy.o: commons.o
unresdump.o: commons.o
unresguessts.o: commons.o key.o modtwoend.o modunres.o
unresnebguessts.o: commons.o key.o modunres.o
unresoptim.o: commons.o modunres.o
unressetdihe.o: modunres.o
unrestransform.o: commons.o modhess.o modunres.o
unrestwist.o: commons.o modunres.o
utils.o: commons.o
vdump.o: key.o modhess.o porfuncs.o
welch.o: modhess.o porfuncs.o
xmylbfgs.o: commons.o key.o modtwoend.o porfuncs.o
Zetterling.o: key.o modhess.o porfuncs.o
Ackland_metals.o: ackland_mendelev_fe.h ackland_mishin_cu.h ackland_sma.h
amb_natinterns.o: commons.o intcommons.o modamber9.o
ambnatintern_extras.o: commons.o intcommons.o key.o modamber9.o modmxatms.o \
	specialfuncts.o
bbrsdm.o: commons.o key.o modcharmm.o
bhinterp.o: commons.o efol.o key.o modamber9.o modcharmm.o porfuncs.o \
	specialfuncts.o
binaryio.o: commons.o key.o porfuncs.o
bisect.o: commons.o key.o modamber9.o modcharmm.o porfuncs.o
bulkmindist.o: key.o
checkedtrig.o: msevb_common.o
chiralhyd.o: commons.o intcommons.o key.o modamber9.o
cleanmemory.o: msevb_common.o
congrad.o: commons.o key.o porfuncs.o
cubsplstring.o: gsdata.o
dbpg.o: commons.o key.o modhess.o
dbptd.o: commons.o key.o modhess.o
drvmsevb.o: commons.o modhess.o msevb_common.o msevb_interfaces.o
efol.o: commons.o key.o modcharmm.o modhess.o porfuncs.o syminf.o zwk.o
fmsevb.o: commons.o msevb_common.o msevb_interfaces.o
gb.o: key.o modhess.o
gbbox.o: commons.o key.o modhess.o
gbd.o: key.o modhess.o
grouprotation.o: commons.o
growstring.o: commons.o cubsplstring.o gsdata.o intcommons.o intcoords.o \
	key.o specialfuncts.o
gsdata.o: commons.o intcommons.o key.o
guesspath.o: key.o modcharmm.o modguess.o modunres.o porfuncs.o
intcommons.o: modamber9.o
intcoords.o: commons.o intcommons.o intcoords.o key.o modamber9.o modcharmm.o \
	porfuncs.o specialfuncts.o
intgradlj.o: commons.o key.o
intlbfgs.o: amhglobals.o commons.o efol.o key.o porfuncs.o
intlbfgslj.o: commons.o efol.o key.o porfuncs.o
ljcapsidmodule.o: commons.o key.o modhess.o pymodule.o
lwotpbox.o: commons.o key.o modhess.o
lwotpgh.o: commons.o key.o modhess.o
main.o: commons.o key.o modamber9.o porfuncs.o
make_conpot.o: commons.o key.o
meccano.o: commons.o key.o modcharmm.o modmec.o modneb.o modtwoend.o \
	modunres.o porfuncs.o
minpermdist.o: intcommons.o intcoords.o key.o modamber9.o modcharmm.o
minpermdistrbcom.o: commons.o key.o
minpermrb.o: commons.o intcommons.o intcoords.o key.o modcharmm.o
msevb.o: commons.o msevb_common.o msevb_interfaces.o
msevb_interfaces.o: msevb_common.o
msstock.o: commons.o key.o modhess.o
mycpu_time.o: key.o porfuncs.o
newcapsid.o: commons.o
newmindist.o: commons.o key.o
newtip.o: commons.o key.o modhess.o
pahagh.o: commons.o key.o modhess.o
patchyd.o: commons.o modhess.o
pyg.o: commons.o key.o modhess.o
rbperm.o: commons.o key.o
rigidb.o: commons.o key.o modhess.o
specialfuncts.o: commons.o
stockghaa.o: commons.o key.o modhess.o
ttm3f.o: ttm3f.o
unguesspath.o: key.o modcharmm.o modguess.o modunres.o porfuncs.o
unmeccano.o: commons.o key.o modcharmm.o modhess.o modmec.o modneb.o \
	modtwoend.o modunres.o porfuncs.o zwk.o
