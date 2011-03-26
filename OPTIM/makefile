
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

ROOT_PATH=$(PWD)/../../
INC_PATH=$(ROOT_PATH)/INCLUDE/

PROGNAME=OPTIM
GENFFILES = porfuncs.f90 header.f90

include $(INC_PATH)/def

# }}}

# Object files {{{

ALLSOURCE := $(wildcard *.f) $(wildcard *.f90) $(wildcard *.F)

##OTUSEDSOURCE: source files not used in present GMIN build {{{
## }}}

# NOTUSEDSOURCE: source files not used in present GMIN build {{{
NOTUSEDSOURCE=
## }}} 

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

# DEPENDENCIES {{{

$(PROG): libnn.a libnc.a $(LLIBS)

libnn.a: modcharmm.o key.o commons.o modunres.o porfuncs.o modmec.o modguess.o efol.o growstring.o gsdata.o intcommons.o intcoords.o amhglobals.o 
libnc.a: libnn.a key.o syminf.o modhess.o amhglobals.o
libamh.a: amhglobals.o

2Dfunc.o:      commons.o modhess.o
Ackland_wrapper.o:      key.o
BLN.o:      key.o modhess.o
Clatmin.o:      commons.o
DF1.o:      key.o modhess.o
GoOptim.o:      key.o
JM.o:      modhess.o porfuncs.o
MB.o:      modhess.o commons.o
PV.o:      commons.o key.o
PachecoC60.o:      modhess.o
SW.o:      key.o modhess.o
SiO2.o:      porfuncs.o modhess.o
TIPmodes.o:      modhess.o
Zetterling.o:      porfuncs.o modhess.o key.o
adm.o:      porfuncs.o commons.o key.o
amber.o:      commons.o modamber.o modamber2.o modhess.o key.o
axdiff.o:      modhess.o
aziz.o:      modhess.o
beig.o:      commons.o key.o modneb.o modtwoend.o porfuncs.o
bfgsts.f:      commons.o key.o vecck.o zwk.o modcharmm.o modhess.o porfuncs.o modtwoend.o 
c10.o:        modhess.o
c60diff.o:      modhess.o
c60p.o:      modhess.o
caardiff.o:      modhess.o
capsid.o:      modhess.o
changep.o:      commons.o key.o modtwoend.o modamber.o modneb.o porfuncs.o
charmmBildc.o:      modamber9.o commons.o key.o
connect.o:      commons.o key.o modtwoend.o modneb.o modcharmm.o modunres.o porfuncs.o porfuncs.o syminf.o modhess.o modcharmm.o modunres.o porfuncs.o porfuncs.o porfuncs.o modneb.o key_neb.o libnn.a modcharmm.o modunres.o efol.o porfuncs.o modcharmm.o modamber9.o 

cpmdlatmin.o:      key.o porfuncs.o
ctest.o:      modhess.o
dctrap.o:      key.o modhess.o
dftb.o:      key.o modhess.o
diis.o:      key.o gdiis.o
dimer.o:       modhess.o commons.o modhess.o commons.o
double.o:      modhess.o
dtrap.o:      modhess.o
dumpit.o:      commons.o key.o
dumpp.o:      commons.o key.o  modcharmm.o porfuncs.o
dzugutov.o:      modhess.o
ectrap.o:      key.o
escp.o:      modhess.o
fd.o:      modhess.o key.o
fetchz.o:      commons.o key.o modtwoend.o syminf.o modneb.o modcharmm.o modamber9.o modunres.o libnc.a modguess.o modmec.o porfuncs.o msevb_common.o binaryio.o amhglobals.o
g2special.o:      modhess.o key.o porfuncs.o
g46merdiff.o:        modhess.o
geopt.o:      commons.o key.o modtwoend.o modhess.o modcharmm.o porfuncs.o binaryio.o efol.o modneb.o amhglobals.o modhess.o modhess.o modhess.o modcharmm.o 
main.o:      commons.o porfuncs.o header.o key.o modamber9.o
gmetry.o:      key.o commons.o
greatcirc.o:      commons.o key.o modtwoend.o modhess.o zwk.o modunres.o modcharmm.o porfuncs.o
h2o.o:      modhess.o
hessout.o:      modhess.o
hybridmin.o:      commons.o key.o vecck.o zwk.o modcharmm.o modhess.o porfuncs.o
inertia.o:      commons.o key.o 
input.o:      porfuncs.o
intbeig.o:      commons.o key.o modneb.o modtwoend.o
intbfgsts.o:      commons.o key.o vecck.o zwk.o modunres.o modhess.o
bfgsts.o:      commons.o key.o  vecck.o zwk.o modcharmm.o modhess.o modtwoend.o porfuncs.o
intsecdiag.o:      commons.o key.o
intxmylbfgs.o:      commons.o key.o modtwoend.o
ions.o:      modhess.o
keyword.o:      commons.o key.o modmec.o modtwoend.o modamber.o modamber9.o modneb.o modcharmm.o modunres.o key_neb.o libnc.a modmec.o modguess.o porfuncs.o gay-berne.o msevb_common.o wc.o binaryio.o gsdata.o cubsplstring.o intcommons.o amhglobals.o modmxatms.o
ljdiff.o:      modhess.o
ljms.o:      modhess.o
ljpbin.o:      key.o modhess.o
ljpdiff.o:      key.o modhess.o
ljpkob.o:      key.o modhess.o
ljpshift.o:      key.o modhess.o
mdiff.o:      modhess.o
mied.o:      modhess.o
mindist.o:      key.o porfuncs.o
morph.o:      commons.o key.o zwk.o modtwoend.o modcharmm.o porfuncs.o
morse.o:      modhess.o
mpdiff.o:      modhess.o key.o
msdiff.o:      modhess.o
mylbfgs.o:      commons.o key.o modtwoend.o modhess.o zwk.o modunres.o modcharmm.o modamber9.o porfuncs.o specialfuncts.o
myorient.o:      key.o
odesd.o:      commons.o key.o modhess.o porfuncs.o
oldneb.o:      commons.o modtwoend.o modneb.o key.o modhess.o porfuncs.o commons.o modcharmm.o
p46merdiff.o:        modhess.o
p4diff.o:      modhess.o
path.o:      commons.o key.o syminf.o modcharmm.o modunres.o modhess.o efol.o porfuncs.o amhglobals.o
pertable.o:      commons.o key.o
potential.o:      commons.o key.o modhess.o modcharmm.o porfuncs.o modamber9.o commons.o modhess.o
projh.o:      modhess.o
rad.o:      commons.o key.o
repelsp.o:      commons.o key.o
rotd.o:      modhess.o
scdiff.o:      modhess.o
secdiag.o:      commons.o key.o modcharmm.o porfuncs.o
shifth.o:      key.o modhess.o
sortxyz.o:      commons.o
stock.o:      key.o modhess.o commons.o vecck.o modhess.o
symmetry.o:      commons.o key.o syminf.o
tcheck.o:      commons.o
thomson.o:      modhess.o commons.o key.o 
tight.o:      key.o modhess.o
tip.o:      commons.o key.o modhess.o
tosifumi.o:      porfuncs.o modhess.o
twoend.o:      porfuncs.o commons.o key.o modtwoend.o zwk.o modcharmm.o modhess.o
unrescalcdihe.f:      commons.o modunres.o commons.o modunres.o modunres.o commons.o
unresconnectsections.o:      commons.o modunres.o commons.o key.o modtwoend.o modunres.o commons.o modtwoend.o modunres.o
unresdummy.o:      commons.o
unresdump.o:      commons.o
unresguessts.o:      commons.o key.o modtwoend.o modunres.o commons.o modtwoend.o modunres.o
unresnebguessts.o:      commons.o key.o modunres.o key_neb.o commons.o modunres.o
unresoptim.o:      commons.o modunres.o commons.o modunres.o commons.o modunres.o commons.o modunres.o
unressetdihe.o:      modunres.o
unrestransform.o:      modunres.o commons.o modhess.o modunres.o modunres.o
unrestwist.o:      commons.o modunres.o commons.o modunres.o commons.o modunres.o
utils.o:      commons.o
vdump.o:      key.o modhess.o porfuncs.o
welch.o:      porfuncs.o modhess.o
xmylbfgs.o:      commons.o key.o modtwoend.o porfuncs.o
amb_natinterns.o:      commons.o modamber9.o intcommons.o 
ambnatintern_extras.o:      modamber9.o intcommons.o modmxatms.o specialfuncts.o commons.o key.o 
bbrsdm.o:      commons.o key.o modcharmm.o
bhinterp.o: key.o modcharmm.o modamber9.o porfuncs.o commons.o efol.o specialfuncts.o
binaryio.o:          key.o commons.o porfuncs.o commons.o
bisect.o: key.o modcharmm.o modamber9.o porfuncs.o commons.o
checkedtrig.o:         msevb_common.o
chiralhyd.o:      commons.o modamber9.o intcommons.o commons.o modamber9.o key.o 
cleanmemory.o:  msevb_common.o
cubsplstring.o:  gsdata.o
dbpg.o:      modhess.o commons.o key.o commons.o 
dbptd.o:     modhess.o commons.o key.o commons.o 
drvmsevb.o:	commons.o msevb_common.o msevb_interfaces.o msevb_interfaces.o modhess.o
efol.o:      commons.o key.o syminf.o modhess.o zwk.o porfuncs.o modcharmm.o
fmsevb.o:	commons.o msevb_common.o msevb_interfaces.o
gay-berne.o: key.o commons.o commons.o modhess.o
gb.o:      modhess.o key.o
gbbox.o:      modhess.o commons.o key.o
gbd.o:      modhess.o key.o
growstring.o:  commons.o key.o cubsplstring.o gsdata.o specialfuncts.o intcommons.o intcoords.o
gsdata.o:  commons.o intcommons.o key.o
guesspath.o: key.o modcharmm.o modguess.o modunres.o porfuncs.o
intcommons.o:  modamber9.o
intcoords.o:  modcharmm.o commons.o intcommons.o key.o specialfuncts.o modamber9.o porfuncs.o
lwotpbox.o: modhess.o commons.o key.o 
lwotpgh.o: modhess.o commons.o key.o 
meccano.o: commons.o porfuncs.o key.o modneb.o modtwoend.o modmec.o key_neb.o commons.o modunres.o modcharmm.o 
meccano.o: key.o
minpermdist.o: key.o modcharmm.o modamber9.o intcommons.o intcoords.o commons.o
minpermrb.o: key.o modcharmm.o intcommons.o commons.o 
msevb.o:     commons.o msevb_common.o msevb_interfaces.o
msevb_interfaces.o:       msevb_common.o
msstock.o:     modhess.o commons.o key.o
mycpu_time.o: porfuncs.o key.o
newcapsid.o:      commons.o 
newmindist.o: commons.o key.o
newtip.o:      commons.o key.o modhess.o
pahagh.o: modhess.o commons.o key.o
patchyd.o: modhess.o commons.o
pyg.o:    modhess.o commons.o key.o
rbperm.o: commons.o key.o
minpermdistrbcom.o: commons.o key.o
rigidb.o: key.o modhess.o commons.o
specialfuncts.o:    commons.o
stockghaa.o: modhess.o commons.o key.o 
unguesspath.o: key.o modcharmm.o modguess.o modunres.o porfuncs.o
unmeccano.o: commons.o key.o modneb.o modunres.o modmec.o key_neb.o commons.o modtwoend.o zwk.o modcharmm.o modneb.o modmec.o porfuncs.o modhess.o
OPTIM.f: OPTIM.F
OPTIM.o: commons.o key.o modtwoend.o modhess.o vecck.o zwk.o modcharmm.o modunres.o libnn.a libnc.a modguess.o modmec.o porfuncs.o charutils.o intcommons.o intcoords.o 
bulkmindist.o:  key.o
fetchz.o: SD.o 
#potential.o: SD.o finite_differences.o libbowman.a
potential.o: SD.o finite_differences.o ttm3f.o VRTMCY5f.o H2O.pjt2.o
make_conpot.o: key.o
grouprotation.o: commons.o
intlbfgs.o: key.o porfuncs.o commons.o intcommons.o amhglobals.o efol.o
intlbfgslj.o: key.o porfuncs.o commons.o intcommons.o efol.o
congrad.o: key.o commons.o 
intgradlj.o: key.o commons.o 
EYtrap.o: key.o
#}}}
