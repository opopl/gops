
WRITE(LFH,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
WRITE(LFH,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'

      DO ISTEP=NDONE+1,NSTEPS 
         CALL FLUSH(LFH)
         MCTEMP = TEMP

         SAVECOORDS=COORDS
 
         DO IA=1,NATOMS
                CALL GET_RND(RND,3,-1.0D0,1.0D0)
                COORDS(IA,1:3)=COORDS(IA,1:3)+STEP*RND(1:3)
         ENDDO
              
         NQ=NQ+1

CALL QUENCH(SCREENC,ITERATIONS,TIME,QDONE)  

WRITE(LFH,111)  'Qu ',          NQ,	&
                ' E=',          POTEL,	&
                ' steps=',      ITERATIONS,	&
                ' RMS=',        RMS,	&
                ' Markov E=',   EPREV,	&
                ' t=',TIME-TSTART	

!
include mc.bh.pairdist.inc.f90
include mc.bh.trackdata.inc.f90

            ELSE
               ATEST=.TRUE. 
               IF (QDONE.EQ.0) THEN
                  ATEST=.FALSE.
               ELSEIF (ATEST) THEN
                  CALL TRANSITION(POTEL,EPREV,ATEST,RANDOM,MCTEMP)
               ENDIF
include mc.bh.checkmarkov.inc.f90

! Accept or reject step
include mc.bh.accrej.inc.f90
! Check the acceptance ratio.
            IF ( (MOD(ISTEP,NACCEPT).EQ.0) CALL ACCREJ(NSUCCESS,NFAIL,NSUCCESST,NFAILT)

            TEMP=TEMP*SCALEFAC
            IF (HIT) GOTO 37
            IF (DUMPINT.GT.0) THEN
               IF (MOD(ISTEP,DUMPINT).EQ.0) THEN
                  CALL DUMPSTATE(ISTEP,EBEST,BESTCOORDS,JBEST)
               ENDIF
            ENDIF
           IF (NQ.GT.NSTEPS) GOTO 37

         ENDDO
 
         CALL FLUSH(LFH)
      ENDDO

