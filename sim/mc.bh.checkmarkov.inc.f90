!
!  Sanity check to make sure the Markov energy agrees with COORDSO. 
!  Stop if not true.
!
               IF (DEBUG.OR.CHECKMARKOVT) THEN
                  CALL POTENTIAL(COORDSO(:),GRAD,OPOTEL,.FALSE.,.FALSE.)
                  IF (ABS(OPOTEL-EPREV).GT.ECONV) THEN
                     IF (EVAP) THEN
                        WRITE(LFH,22) 
                        'mc> WARNING - energy for saved coordinates ',OPOTEL,&
                        ' differs from Markov energy ',
                        EPREV,' because an atom moved outside the container'
                     ELSE
                        WRITE(LFH,'(2(A,G20.10))') 'mc> ERROR - energy for coordinates in COORDSO=',OPOTEL,
     &                                                 ' but Markov energy=',EPREV
                        STOP
                     ENDIF
                  ENDIF
               ENDIF

