# Project dir: /home/gates/gops/O
# Program name: O
# /home/gates/gops/O/deps.mk
# Fortran dependency file
# Created: 15:49:4, Sat Oct 1, 2011
Clatmin.o: commons.o
stock.o: commons.o key.o modhess.o vecck.o
intsecdiag.o: commons.o key.o
fd.o: key.o modhess.o
gmetry.o: commons.o key.o
ions.o: modhess.o
potential.o: SD.o VRTMCY5f.o commons.o finite_differences.o key.o modamber9.o \
	modcharmm.o modhess.o porfuncs.o
secdiag.o: commons.o key.o modcharmm.o porfuncs.o
hessout.o: modhess.o
ttm3f.o: ttm3f.o
SiO2.o: modhess.o porfuncs.o
myorient.o: key.o
congrad.o: commons.o key.o porfuncs.o
mied.o: modhess.o
binaryio.o: commons.o key.o porfuncs.o
msevb_interfaces.o: msevb_common.o
vdump.o: key.o modhess.o porfuncs.o
scdiff.o: modhess.o
patchyd.o: commons.o modhess.o
cleanmemory.o: msevb_common.o
bhinterp.o: commons.o efol.o key.o modamber9.o modcharmm.o porfuncs.o \
	specialfuncts.o
unresguessts.o: commons.o key.o modtwoend.o modunres.o vars.o
newtip.o: commons.o key.o modhess.o
intbeig.o: commons.o key.o modneb.o modtwoend.o
xmylbfgs.o: commons.o key.o modtwoend.o porfuncs.o
intlbfgs.o: amhglobals.o commons.o efol.o key.o porfuncs.o
ljms.o: modhess.o
mindist.o: key.o porfuncs.o
pahagh.o: commons.o key.o modhess.o
adm.o: commons.o key.o porfuncs.o
intxmylbfgs.o: commons.o key.o modtwoend.o porfuncs.o
msevb.o: commons.o msevb_common.o msevb_interfaces.o
intgradlj.o: commons.o key.o
drvmsevb.o: commons.o modhess.o msevb_common.o msevb_interfaces.o
mylbfgs.o: commons.o key.o modamber9.o modcharmm.o modhess.o modtwoend.o \
	modunres.o porfuncs.o specialfuncts.o zwk.o
morph.o: commons.o key.o modcharmm.o modtwoend.o porfuncs.o zwk.o
fmsevb.o: commons.o msevb_common.o msevb_interfaces.o
tip.o: commons.o key.o modhess.o
thomson.o: commons.o key.o modhess.o
path.o: amhglobals.o commons.o efol.o key.o modcharmm.o modhess.o modunres.o \
	porfuncs.o syminf.o
newmindist.o: commons.o key.o
ljpbin.o: key.o modhess.o
minpermrb.o: commons.o intcommons.o intcoords.o key.o modcharmm.o
keyword.o: amhglobals.o binaryio.o commons.o cubsplstring.o gsdata.o \
	intcommons.o key.o libnc.a libnn.a modamber.o modamber9.o modcharmm.o \
	modguess.o modmec.o modmxatms.o modneb.o modtwoend.o modunres.o \
	msevb_common.o porfuncs.o pymodule.o wc.o
dbptd.o: commons.o key.o modhess.o
c60diff.o: modhess.o
axdiff.o: modhess.o
ctest.o: modhess.o
input.o: porfuncs.o
unresnebguessts.o: commons.o key.o libnn.a modunres.o
odesd.o: commons.o key.o modhess.o porfuncs.o
unresconnectsections.o: commons.o key.o modtwoend.o modunres.o vars.o
gsdata.o: commons.o intcommons.o key.o
meccano.o: commons.o key.o libnn.a modcharmm.o modmec.o modneb.o modtwoend.o \
	modunres.o porfuncs.o
mdiff.o: modhess.o
shifth.o: key.o modhess.o
tight.o: key.o modhess.o
Ackland_wrapper.o: key.o
twoend.o: commons.o key.o modcharmm.o modhess.o modtwoend.o porfuncs.o zwk.o
morse.o: modhess.o
amber.o: commons.o key.o modamber.o modamber2.o modhess.o
p4diff.o: modhess.o
OPTIM.o: commons.o intcommons.o intcoords.o key.o libnc.a libnn.a modcharmm.o \
	modguess.o modhess.o modmec.o modneb.o modtwoend.o modunres.o \
	porfuncs.o vecck.o zwk.o
EYtrap.o: key.o modhess.o
ljcapsidmodule.o: commons.o key.o modhess.o pymodule.o
BLN.o: key.o modhess.o
mycpu_time.o: key.o porfuncs.o
gbbox.o: commons.o key.o modhess.o
rca.o: commons.o dv.o f.o porfuncs.o v.o
inertia.o: commons.o key.o
dctrap.o: key.o modhess.o
DF1.o: key.o modhess.o
dzugutov.o: modhess.o
stockghaa.o: commons.o key.o modhess.o
intbfgsts.o: commons.o key.o modhess.o modunres.o porfuncs.o vecck.o zwk.o
c10.o: modhess.o
changep.o: commons.o key.o modamber.o modneb.o modtwoend.o porfuncs.o
mpdiff.o: key.o modhess.o
ectrap.o: key.o
GoOptim.o: key.o
minpermdistrbcom.o: commons.o key.o
repelsp.o: commons.o key.o
msdiff.o: modhess.o
connect.o: commons.o efol.o key.o libnn.a modamber9.o modcharmm.o modhess.o \
	modneb.o modtwoend.o modunres.o porfuncs.o syminf.o
caardiff.o: modhess.o
unresdump.o: commons.o
dbpg.o: commons.o key.o modhess.o
MB.o: commons.o modhess.o
cpmdlatmin.o: key.o porfuncs.o
lwotpbox.o: commons.o key.o modhess.o
diis.o: gdiis.o key.o
ljpkob.o: key.o modhess.o
PachecoC60.o: modhess.o
grouprotation.o: commons.o
tcheck.o: commons.o
cubsplstring.o: gsdata.o
gb.o: key.o modhess.o
unrestransform.o: commons.o modhess.o modunres.o
pertable.o: commons.o key.o
checkedtrig.o: msevb_common.o
minpermdist.o: intcommons.o intcoords.o key.o modamber9.o modcharmm.o
double.o: modhess.o
rotd.o: modhess.o
dtrap.o: modhess.o
newcapsid.o: commons.o
dumpit.o: commons.o key.o
lwotpgh.o: commons.o key.o modhess.o
dimer.o: commons.o modhess.o
charmmBildc.o: commons.o modamber9.o
unguesspath.o: key.o modcharmm.o modguess.o modunres.o porfuncs.o
sortxyz.o: commons.o
amb_natinterns.o: commons.o intcommons.o modamber9.o
unmeccano.o: commons.o key.o libnn.a modcharmm.o modhess.o modmec.o modneb.o \
	modtwoend.o modunres.o porfuncs.o zwk.o
growstring.o: commons.o cubsplstring.o gsdata.o intcommons.o intcoords.o \
	key.o specialfuncts.o
PV.o: commons.o key.o
h2o.o: modhess.o
2Dfunc.o: commons.o modhess.o
dftb.o: key.o modhess.o
bulkmindist.o: key.o
bisect.o: commons.o key.o modamber9.o modcharmm.o porfuncs.o
symmetry.o: commons.o key.o syminf.o
gbd.o: key.o modhess.o
welch.o: modhess.o porfuncs.o
unrescalcdihe.o: commons.o modunres.o
unresoptim.o: commons.o modunres.o
ljpdiff.o: key.o modhess.o
rigidb.o: commons.o key.o modhess.o
rbperm.o: commons.o key.o
Zetterling.o: key.o modhess.o porfuncs.o
aziz.o: modhess.o
beig.o: commons.o key.o modneb.o modtwoend.o porfuncs.o
unressetdihe.o: modunres.o
ljdiff.o: modhess.o
geopt.o: amhglobals.o binaryio.o commons.o efol.o key.o modcharmm.o modhess.o \
	modneb.o modtwoend.o porfuncs.o
msstock.o: commons.o key.o modhess.o
hybridmin.o: commons.o key.o modcharmm.o modhess.o porfuncs.o vecck.o zwk.o
intcommons.o: modamber9.o
utils.o: commons.o
ambnatintern_extras.o: commons.o intcommons.o key.o modamber9.o modmxatms.o \
	specialfuncts.o
rad.o: commons.o key.o
SW.o: key.o modhess.o
oldneb.o: commons.o key.o modcharmm.o modhess.o modneb.o modtwoend.o \
	porfuncs.o
specialfuncts.o: commons.o
g46merdiff.o: modhess.o
greatcirc.o: commons.o key.o modcharmm.o modhess.o modtwoend.o modunres.o \
	porfuncs.o zwk.o
efol.o: commons.o key.o modcharmm.o modhess.o porfuncs.o syminf.o zwk.o
TIPmodes.o: modhess.o
bbrsdm.o: commons.o key.o modcharmm.o
make_conpot.o: commons.o key.o
intcoords.o: commons.o intcommons.o intcoords.o key.o modamber9.o modcharmm.o \
	porfuncs.o specialfuncts.o
unrestwist.o: commons.o modunres.o
chiralhyd.o: commons.o intcommons.o key.o modamber9.o
fetchz.o: SD.o amhglobals.o binaryio.o commons.o key.o libnc.a modamber9.o \
	modcharmm.o modguess.o modmec.o modneb.o modtwoend.o modunres.o \
	msevb_common.o porfuncs.o syminf.o
p46merdiff.o: modhess.o
tosifumi.o: modhess.o porfuncs.o
escp.o: modhess.o
capsid.o: modhess.o
JM.o: modhess.o porfuncs.o
projh.o: modhess.o
dumpp.o: commons.o key.o modcharmm.o porfuncs.o
bfgsts.o: commons.o key.o modcharmm.o modhess.o modtwoend.o porfuncs.o \
	vecck.o zwk.o
ljpshift.o: key.o modhess.o
c60p.o: modhess.o
intlbfgslj.o: commons.o efol.o key.o porfuncs.o
g2special.o: key.o modhess.o porfuncs.o
main.o: commons.o key.o modamber9.o porfuncs.o
guesspath.o: key.o modcharmm.o modguess.o modunres.o porfuncs.o
pyg.o: commons.o key.o modhess.o
