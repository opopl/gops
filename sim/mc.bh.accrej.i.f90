! Accept or reject step. If the quench did not converge then allow a
! potential move, but count it as a rejection in terms of NSUCCESS and
! NFAIL. This way we will accept a lower minimum if found, 
! but the steps won't become so big.
               IF (ATEST) THEN
                  IF ((ISTEP-JACCPREV.GT.NRELAX).AND.ABS(DQE).GT.EDIFF) THEN
                     JACCPREV=ISTEP
                  ENDIF
                  IF (QDONE.EQ.0) THEN
                     NSUCCESS=NSUCCESS+1
                  ELSE
                     NFAIL=NFAIL+1
                  ENDIF
                  QEPPREV=QEPREV
                  QEPREV=QE
               ELSE
                  NFAIL=NFAIL+1
                  !CALL MYRESET(NATOMS,NSEED)
               ENDIF

