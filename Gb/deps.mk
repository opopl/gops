# Project dir: /home/op/gops/Gb
# Not used dir: old
# Not used dir: c
# Not used dir: n
# Not used dir: nu
# /home/op/gops/Gb/deps.mk
# Fortran dependency file
# Created: 3:4:10, Wed Sep 21, 2011
finalq.o: commons.o qmod.o
g46merdiff.o: modhess.o
commons.o: countatoms.o
mcruns.o: commons.o mc.o
centre.o: commons.o
output.o: commons.o file_manager.o
io1.o: commons.o porfuncs.o qmod.o
main.o: commons.o countatoms.o mc.o porfuncs.o qmod.o
sort2.o: qmod.o
mc.o: commons.o porfuncs.o qmod.o
takestep.o: commons.o
centrecom.o: commons.o
GMINdump.o: commons.o output.o porfuncs.o qmod.o
mylbfgs.o: commons.o porfuncs.o
finalio.o: commons.o qmod.o
p46merdiff.o: modhess.o
hsmove.o: commons.o
rca.o: commons.o dv.o porfuncs.o
potential.o: commons.o porfuncs.o qmod.o
mycpu_time.o: commons.o
quench.o: commons.o porfuncs.o qmod.o
BLN.o: commons.o
