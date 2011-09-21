# Project dir: /home/op/gops/Gb
# Not used dir: old
# Not used dir: c
# Not used dir: n
# Not used dir: nu
# /home/op/gops/Gb/deps.mk
# Fortran dependency file
# Created: 1:55:8, Wed Sep 21, 2011
finalq.o: qmod.o v.o
mcruns.o: v.o
io1.o: porfuncs.o qmod.o v.o
main.o: countatoms.o porfuncs.o qmod.o v.o
mc.o: porfuncs.o qmod.o v.o
takestep.o: v.o
mylbfgs.o: porfuncs.o v.o
finalio.o: qmod.o v.o
rca.o: dv.o porfuncs.o v.o
potential.o: porfuncs.o qmod.o v.o
v.o: countatoms.o
quench.o: porfuncs.o qmod.o v.o
