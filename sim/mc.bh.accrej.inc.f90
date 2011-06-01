! Accept or reject step. If the quench did not converge then allow a
! potential move, but count it as a rejection in terms of NSUCCESS and
! NFAIL. This way we will accept a lower minimum if found, 
! but the steps won't become so big.
               IF (ATEST) THEN
                  IF (DEBUG) WRITE(LFH,34) JP,RANDOM,POTEL,EPREV,NSUCCESS,NFAIL
                  IF ((ISTEP-JACCPREV.GT.NRELAX).AND.ABS(POTEL-EPREV).GT.ECONV) THEN
                     JACCPREV=ISTEP
                  ENDIF
                  IF (QDONE.EQ.0) THEN
                     NSUCCESS=NSUCCESS+1
                  ELSE
                     NFAIL=NFAIL+1
                  ENDIF
                  EPPREV=EPREV
                  EPREV=POTEL
                  COORDSO=COORDS
                  VATO=VAT
               ELSE
                  NFAIL=NFAIL+1
                  CALL MYRESET(NATOMS,NSEED)
                  IF (DEBUG) THEN
                     WRITE(LFH,36) JP,RANDOM,POTEL,EPREV,NSUCCESS,NFAIL
36                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' REJ')
                  ENDIF
               ENDIF
            ENDIF

