
WRITE(LFH,'(A)') 'Calculating initial energy'

FQFLAG=.FALSE.

CALL QUENCH(SCREENC,ITERATIONS,TIME,QDONE)

include write.111.inc.f90

EPREV=E

