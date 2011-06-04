
S(6)=S(1)**6 ; S(12)=S(6)**2
GNB= 0.0D0; GB= 0.0D0 ; GBA= 0.0D0 ; GTA= 0.0D0 ; GRAD= 0.0D0 

include bln.grad.nbond.inc.f90  ! GNB         non-bonded 
include bln.grad.bond.inc.f90   ! GB          bonded  
include bln.grad.bangle.inc.f90 ! GBA         bond angles
include bln.grad.tangle.inc.f90 ! GTA         torsional angles

G=GNB+GB+GBA+GTA 

