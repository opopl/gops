
S(6)=S(1)**6 ; S(12)=S(6)**2
GNB= 0.0D0; GB= 0.0D0 ; GBA= 0.0D0 ; GTA= 0.0D0 ; GRAD= 0.0D0 

include bln.grad.nbond.inc.f90
include bln.grad.bond.inc.f90
include bln.grad.bangle.inc.f90
include bln.grad.tangle.inc.f90

GRAD=GNB+GB+GBA+GTA 

