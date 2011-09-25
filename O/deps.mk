# Project dir: /home/op226/gops/O
# Program name: O
# /home/op226/gops/O/deps.mk
# Fortran dependency file
# Created: 20:13:53, Sun Sep 25, 2011
2Dfunc.o: commons.o modhess.o
Ackland_wrapper.o: key.o
BLN.o: key.o modhess.o
Clatmin.o: commons.o
DF1.o: key.o modhess.o
EYtrap.o: key.o modhess.o
GoOptim.o: key.o
JM.o: modhess.o porfuncs.o
MB.o: commons.o modhess.o
OPTIM.o: commons.o intcommons.o intcoords.o key.o libnc.a libnn.a modcharmm.o \
	modguess.o modhess.o modmec.o modneb.o modtwoend.o modunres.o \
	porfuncs.o vecck.o zwk.o
PV.o: commons.o key.o
PachecoC60.o: modhess.o
SW.o: key.o modhess.o
SiO2.o: modhess.o porfuncs.o
TIPmodes.o: modhess.o
Zetterling.o: key.o modhess.o porfuncs.o
adm.o: commons.o key.o porfuncs.o
amb_natinterns.o: commons.o intcommons.o modamber9.o
amber.o: commons.o key.o modamber.o modamber2.o modhess.o
ambnatintern_extras.o: commons.o intcommons.o key.o modamber9.o modmxatms.o \
	specialfuncts.o
axdiff.o: modhess.o
aziz.o: modhess.o
bbrsdm.o: commons.o key.o modcharmm.o
beig.o: commons.o key.o modneb.o modtwoend.o porfuncs.o
bfgsts.o: commons.o key.o modcharmm.o modhess.o modtwoend.o porfuncs.o \
	vecck.o zwk.o
bhinterp.o: commons.o efol.o key.o modamber9.o modcharmm.o porfuncs.o \
	specialfuncts.o
binaryio.o: commons.o key.o porfuncs.o
bisect.o: commons.o key.o modamber9.o modcharmm.o porfuncs.o
bulkmindist.o: key.o
c10.o: modhess.o
c60diff.o: modhess.o
c60p.o: modhess.o
caardiff.o: modhess.o
capsid.o: modhess.o
changep.o: commons.o key.o modamber.o modneb.o modtwoend.o porfuncs.o
charmmBildc.o: commons.o modamber9.o
checkedtrig.o: msevb_common.o
chiralhyd.o: commons.o intcommons.o key.o modamber9.o
cleanmemory.o: msevb_common.o
congrad.o: commons.o key.o porfuncs.o
connect.o: commons.o efol.o key.o libnn.a modamber9.o modcharmm.o modhess.o \
	modneb.o modtwoend.o modunres.o porfuncs.o syminf.o
cpmdlatmin.o: key.o porfuncs.o
ctest.o: modhess.o
cubsplstring.o: gsdata.o
dbpg.o: commons.o key.o modhess.o
dbptd.o: commons.o key.o modhess.o
dctrap.o: key.o modhess.o
dftb.o: key.o modhess.o
diis.o: gdiis.o key.o
dimer.o: commons.o modhess.o
double.o: modhess.o
drvmsevb.o: commons.o modhess.o msevb_common.o msevb_interfaces.o
dtrap.o: modhess.o
dumpit.o: commons.o key.o
dumpp.o: commons.o key.o modcharmm.o porfuncs.o
dzugutov.o: modhess.o
ectrap.o: key.o
efol.o: commons.o key.o modcharmm.o modhess.o porfuncs.o syminf.o zwk.o
escp.o: modhess.o
fd.o: key.o modhess.o
fetchz.o: SD.o amhglobals.o binaryio.o commons.o key.o libnc.a modamber9.o \
	modcharmm.o modguess.o modmec.o modneb.o modtwoend.o modunres.o \
	msevb_common.o porfuncs.o syminf.o
fmsevb.o: commons.o msevb_common.o msevb_interfaces.o
g2special.o: key.o modhess.o porfuncs.o
g46merdiff.o: modhess.o
gb.o: key.o modhess.o
gbbox.o: commons.o key.o modhess.o
gbd.o: key.o modhess.o
geopt.o: amhglobals.o binaryio.o commons.o efol.o key.o modcharmm.o modhess.o \
	modneb.o modtwoend.o porfuncs.o
gmetry.o: commons.o key.o
greatcirc.o: commons.o key.o modcharmm.o modhess.o modtwoend.o modunres.o \
	porfuncs.o zwk.o
grouprotation.o: commons.o
growstring.o: commons.o cubsplstring.o gsdata.o intcommons.o intcoords.o \
	key.o specialfuncts.o
gsdata.o: commons.o intcommons.o key.o
guesspath.o: key.o modcharmm.o modguess.o modunres.o porfuncs.o
h2o.o: modhess.o
hessout.o: modhess.o
hybridmin.o: commons.o key.o modcharmm.o modhess.o porfuncs.o vecck.o zwk.o
inertia.o: commons.o key.o
input.o: porfuncs.o
intbeig.o: commons.o key.o modneb.o modtwoend.o
intbfgsts.o: commons.o key.o modhess.o modunres.o porfuncs.o vecck.o zwk.o
intcommons.o: modamber9.o
intcoords.o: commons.o intcommons.o intcoords.o key.o modamber9.o modcharmm.o \
	porfuncs.o specialfuncts.o
intgradlj.o: commons.o key.o
intlbfgs.o: amhglobals.o commons.o efol.o key.o porfuncs.o
intlbfgslj.o: commons.o efol.o key.o porfuncs.o
intsecdiag.o: commons.o key.o
intxmylbfgs.o: commons.o key.o modtwoend.o porfuncs.o
ions.o: modhess.o
keyword.o: amhglobals.o binaryio.o commons.o cubsplstring.o gsdata.o \
	intcommons.o key.o libnc.a libnn.a modamber.o modamber9.o modcharmm.o \
	modguess.o modmec.o modmxatms.o modneb.o modtwoend.o modunres.o \
	msevb_common.o porfuncs.o pymodule.o wc.o
ljcapsidmodule.o: commons.o key.o modhess.o pymodule.o
ljdiff.o: modhess.o
ljms.o: modhess.o
ljpbin.o: key.o modhess.o
ljpdiff.o: key.o modhess.o
ljpkob.o: key.o modhess.o
ljpshift.o: key.o modhess.o
lwotpbox.o: commons.o key.o modhess.o
lwotpgh.o: commons.o key.o modhess.o
main.o: commons.o key.o modamber9.o porfuncs.o
make_conpot.o: commons.o key.o
mdiff.o: modhess.o
meccano.o: commons.o key.o libnn.a modcharmm.o modmec.o modneb.o modtwoend.o \
	modunres.o porfuncs.o
mied.o: modhess.o
mindist.o: key.o porfuncs.o
minpermdist.o: intcommons.o intcoords.o key.o modamber9.o modcharmm.o
minpermdistrbcom.o: commons.o key.o
minpermrb.o: commons.o intcommons.o intcoords.o key.o modcharmm.o
rad.o: commons.o key.o
morph.o: commons.o key.o modcharmm.o modtwoend.o porfuncs.o zwk.o
morse.o: modhess.o
mpdiff.o: key.o modhess.o
msdiff.o: modhess.o
msevb.o: commons.o msevb_common.o msevb_interfaces.o
msevb_interfaces.o: msevb_common.o
msstock.o: commons.o key.o modhess.o
mycpu_time.o: key.o porfuncs.o
mylbfgs.o: commons.o key.o modamber9.o modcharmm.o modhess.o modtwoend.o \
	modunres.o porfuncs.o specialfuncts.o zwk.o
myorient.o: key.o
newcapsid.o: commons.o
newmindist.o: commons.o key.o
newtip.o: commons.o key.o modhess.o
odesd.o: commons.o key.o modhess.o porfuncs.o
oldneb.o: commons.o key.o modcharmm.o modhess.o modneb.o modtwoend.o \
	porfuncs.o
p46merdiff.o: modhess.o
p4diff.o: modhess.o
pahagh.o: commons.o key.o modhess.o
patchyd.o: commons.o modhess.o
path.o: amhglobals.o commons.o efol.o key.o modcharmm.o modhess.o modunres.o \
	porfuncs.o syminf.o
pertable.o: commons.o key.o
potential.o: SD.o VRTMCY5f.o commons.o finite_differences.o key.o modamber9.o \
	modcharmm.o modhess.o porfuncs.o
projh.o: modhess.o
pyg.o: commons.o key.o modhess.o
rbperm.o: commons.o key.o
repelsp.o: commons.o key.o
rigidb.o: commons.o key.o modhess.o
rotd.o: modhess.o
scdiff.o: modhess.o
secdiag.o: commons.o key.o modcharmm.o porfuncs.o
shifth.o: key.o modhess.o
sortxyz.o: commons.o
specialfuncts.o: commons.o
stock.o: commons.o key.o modhess.o vecck.o
stockghaa.o: commons.o key.o modhess.o
symmetry.o: commons.o key.o syminf.o
tcheck.o: commons.o
thomson.o: commons.o key.o modhess.o
tight.o: key.o modhess.o
tip.o: commons.o key.o modhess.o
tosifumi.o: modhess.o porfuncs.o
ttm3f.o: ttm3f.o
twoend.o: commons.o key.o modcharmm.o modhess.o modtwoend.o porfuncs.o zwk.o
unguesspath.o: key.o modcharmm.o modguess.o modunres.o porfuncs.o
unmeccano.o: commons.o key.o libnn.a modcharmm.o modhess.o modmec.o modneb.o \
	modtwoend.o modunres.o porfuncs.o zwk.o
unrescalcdihe.o: commons.o modunres.o
unresconnectsections.o: commons.o key.o modtwoend.o modunres.o vars.o
unresdump.o: commons.o
unresguessts.o: commons.o key.o modtwoend.o modunres.o vars.o
unresnebguessts.o: commons.o key.o libnn.a modunres.o
unresoptim.o: commons.o modunres.o
unressetdihe.o: modunres.o
unrestransform.o: commons.o modhess.o modunres.o
unrestwist.o: commons.o modunres.o
utils.o: commons.o
vdump.o: key.o modhess.o porfuncs.o
welch.o: modhess.o porfuncs.o
xmylbfgs.o: commons.o key.o modtwoend.o porfuncs.o
rca.o: commons.o dv.o porfuncs.o v.o
