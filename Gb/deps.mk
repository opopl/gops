# Project dir: /home/op226/gops/Gb
# Not used dir: old
# Not used dir: c
# Not used dir: n
# Not used dir: nu
# /home/op226/gops/Gb/deps.mk
# Fortran dependency file
# Created: 14:12:1, Wed Sep 21, 2011
output.o: commons.o file_manager.o
GMINdump.o: commons.o output.o porfuncs.o qmod.o
finalq.o: commons.o qmod.o
strings.o: precmod.o
mcruns.o: commons.o mc.o
countatoms.o: commons.o
rca.o: commons.o dv.o func.o porfuncs.o
centre.o: commons.o
quench.o: commons.o porfuncs.o qmod.o
io1.o: commons.o porfuncs.o qmod.o
kw.o: commons.o func.o strings.o
func.o: commons.o porfuncs.o
sort2.o: qmod.o
mc.o: commons.o porfuncs.o qmod.o
centrecom.o: commons.o
BLN.o: commons.o
p46merdiff.o: modhess.o
mylbfgs.o: commons.o porfuncs.o
finalio.o: commons.o qmod.o
mycpu_time.o: commons.o
potential.o: commons.o porfuncs.o qmod.o
g46merdiff.o: modhess.o
main.o: commons.o func.o kw.o mc.o porfuncs.o qmod.o
takestep.o: commons.o
