# Project dir: /home/op226/gops/B
# Program name: B
# /home/op226/gops/B/deps.mk
# Fortran dependency file
# Created: 23:48:49, Fri Sep 23, 2011
BLN.o: commons.o
centre.o: commons.o
centrecom.o: commons.o
convert.o: porfuncs.o
finalq.o: commons.o v.o
g46merdiff.o: v.o
p46merdiff.o: v.o
sort2.o: commons.o v.o
GMINdump.o: commons.o output.o porfuncs.o v.o
finalio.o: commons.o v.o
iom.o: commons.o porfuncs.o v.o
kw.o: commons.o porfuncs.o v.o
main.o: commons.o f.o porfuncs.o v.o
mycpu_time.o: commons.o
mylbfgs.o: commons.o f.o porfuncs.o v.o
myreset.o: commons.o
newmindist.o: commons.o
mc.o: commons.o f.o porfuncs.o v.o
output.o: commons.o file_manager.o
potential.o: bln.o commons.o f.o porfuncs.o v.o
quench.o: commons.o f.o porfuncs.o v.o
saveit.o: commons.o v.o
takestep.o: commons.o
transition.o: commons.o f.o v.o
accrej.o: commons.o v.o
f.o: commons.o v.o
rca.o: commons.o dv.o porfuncs.o
