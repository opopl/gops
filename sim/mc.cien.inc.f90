
WRITE(LFH,'(A)') 'Calculating initial energy'
CALL QUENCH(.FALSE.,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
NQTOT=NQTOT+1
WRITE(LFH,111)  'Qu ',                  NQ,                        	&
                ' E=',                  POTEL,	&
                ' steps=',              ITERATIONS,	&
                ' RMS=',                RMS,	&
                ' Markov E=',           POTEL,	&
                ' t=',                  TIME-TSTART

EPREV=POTEL
EPPREV=0.0D0
EBEST=POTEL
BESTCOORDS=COORDS
JBEST=0
RMIN=POTEL
RCOORDS=COORDS
COORDSO=COORDS
EPSSPHERE=EPSSAVE

