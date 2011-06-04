
WRITE(LFH,'(A)') 'Calculating initial energy'

FQFLAG=.FALSE.

CALL QUENCH(SCREENC,ITERATIONS,TIME,QDONE)

include write.111.inc.f90

!  EPREV saves the previous energy in the Markov chain.
!  EBEST and JBEST record the lowest energy since the last reseeding and the
!  step it was attained at. BESTCOORDS contains the corresponding coordinates.

EPREV=E
EPPREV=ZERO
EBEST=E
JBEST=0


