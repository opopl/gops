# Project dir: /home/op226/gops/Gb
# Not used dir: old
# Not used dir: c
# /home/op226/gops/Gb/deps.mk
# Fortran dependency file
# Created: 0:3:53, Wed Sep 21, 2011
finalq.o: commons.o qmod.o
mcruns.o: commons.o
rca.o: commons.o dv.o porfuncs.o
quench.o: commons.o porfuncs.o qmod.o
io1.o: commons.o porfuncs.o qmod.o
mc1.o: commons.o
mc.o: commons.o porfuncs.o qmod.o
mylbfgs.o: commons.o porfuncs.o
finalio.o: commons.o qmod.o
potential.o: commons.o permu.o porfuncs.o qmod.o
main.o: commons.o countatoms.o porfuncs.o qmod.o
commons.o: countatoms.o
takestep.o: commons.o
