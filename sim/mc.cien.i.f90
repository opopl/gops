
WRITE(LFH,'(A)') 'Calculating initial energy'

FQFLAG=.FALSE.

CALL QUENCH(SCREENC,ITERATIONS,TIME,CONVG)

include "mc.w.111.i.f90"

!  EPREV saves the previous energy in the Markov chain.
!  EBEST and JBEST record the lowest energy since the last reseeding and the
!  step it was attained at. BESTCOORDS contains the corresponding coordinates.

QEPREV=QE
QEPPREV=ZERO
EBEST=QE
JBEST=0


