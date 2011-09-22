# Project dir: /home/op226/gops/Gb
# Program name: Gb
# Not used dir: old
# Not used dir: c
# Not used dir: n
# Not used dir: nu
# /home/op226/gops/Gb/deps.mk
# Fortran dependency file
# Created: 13:46:24, Thu Sep 22, 2011
centre.o: commons.o
potential.o: commons.o porfuncs.o qmod.o
countatoms.o: commons.o
mc.o: commons.o func.o porfuncs.o qmod.o
kw.o: commons.o func.o strings.o
centrecom.o: commons.o
quench.o: commons.o porfuncs.o qmod.o
BLN.o: commons.o
mcruns.o: commons.o mc.o
sort2.o: qmod.o
takestep.o: commons.o
io.o: commons.o porfuncs.o qmod.o
mycpu_time.o: commons.o
output.o: commons.o file_manager.o
main.o: commons.o func.o kw.o mc.o porfuncs.o qmod.o
strings.o: precmod.o
finalio.o: commons.o qmod.o
GMINdump.o: commons.o output.o porfuncs.o qmod.o
finalq.o: commons.o qmod.o
p46merdiff.o: modhess.o
g46merdiff.o: modhess.o
mylbfgs.o: commons.o porfuncs.o
rca.o: commons.o dv.o func.o porfuncs.o
func.o: commons.o porfuncs.o
