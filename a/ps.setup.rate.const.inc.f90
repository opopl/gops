
!  Rate constants.
!op226 {{{
!     KMEAN=0.0D0
      IF (ENSEMBLE.EQ.'T') THEN
        ! {{{
         DO J1=1,NTS
            KPLUS(J1)  = LOG(1.0D0 * HORDERMIN(PLUS(J1))  / (2.0D0 * PI*HORDERTS(J1))) +
     1             (FVIBMIN(PLUS(J1))  - FVIBTS(J1)) / 2.0D0 - (ETS(J1) - EMIN(PLUS(J1)) )/TEMPERATURE
            IF (FRICTIONT) KPLUS(J1)=KPLUS(J1)+LOG(FRICTIONFAC(NEGEIG(J1)))
            KMINUS(J1) = LOG(1.0D0 * HORDERMIN(MINUS(J1)) / (2.0D0 * PI*HORDERTS(J1))) +
     1             (FVIBMIN(MINUS(J1)) - FVIBTS(J1)) / 2.0D0 - (ETS(J1) - EMIN(MINUS(J1)))/TEMPERATURE
            IF (FRICTIONT) KMINUS(J1)=KMINUS(J1)+LOG(FRICTIONFAC(NEGEIG(J1)))
            IF (ZSYM(1:2).EQ.'CA') KPLUS(J1)=KPLUS(J1)+30.66356D0
            IF (ZSYM(1:2).EQ.'CA') KMINUS(J1)=KMINUS(J1)+30.66356D0
            IF (PLUS(J1).EQ.MINUS(J1)) KPLUS(J1)=KPLUS(J1)+LOG(2.0D0)
            IF (PLUS(J1).EQ.MINUS(J1)) KMINUS(J1)=KMINUS(J1)+LOG(2.0D0)
!           KMEAN=KMEAN+KPLUS(J1)+KMINUS(J1)
            IF (DEBUG) WRITE(*,'(A,3I6,5F15.5,3G20.10)') 'setup> J1,PLUS,MINUS,Ets,E+,E-,k+,k-,<k>=',J1,PLUS(J1),MINUS(J1),
     1                                            ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1)
         ENDDO
         ! }}}
      ELSE
        ! {{{
         DO J1=1,NTS
            IF (TOTALE.GT.ETS(J1)) THEN
               KPLUS(J1)  = LOG(1.0D0 * HORDERMIN(PLUS(J1))  / (2*PI*HORDERTS(J1))) +
     1                   (FVIBMIN(PLUS(J1))  - FVIBTS(J1))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J1))/(TOTALE-EMIN(PLUS(J1))))
               KMINUS(J1) = LOG(1.0D0 * HORDERMIN(MINUS(J1)) / (2*PI*HORDERTS(J1))) +
     1                   (FVIBMIN(MINUS(J1)) - FVIBTS(J1))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J1))/(TOTALE-EMIN(MINUS(J1))))
               IF (ZSYM(1:2).EQ.'CA') KPLUS(J1)=KPLUS(J1)+30.66356D0
               IF (ZSYM(1:2).EQ.'CA') KMINUS(J1)=KMINUS(J1)+30.66356D0
               IF (PLUS(J1).EQ.MINUS(J1)) KPLUS(J1)=KPLUS(J1)+LOG(2.0D0)
               IF (PLUS(J1).EQ.MINUS(J1)) KMINUS(J1)=KMINUS(J1)+LOG(2.0D0)
!              KMEAN=KMEAN+KPLUS(J1)+KMINUS(J1)
            ELSE
               KPLUS(J1)=-1.0D250
               KMINUS(J1)=-1.0D250
            ENDIF
         ENDDO
         ! }}}
      ENDIF
      !op226 }}}
!     IF (NTS.GT.0) KMEAN=KMEAN/(2.0D0*NTS)
!     PRINT '(A,G20.10)', 'setup> Mean log rate constant=', KMEAN
!
!  Sums of rates out of the intermediate minima
!
!op226 {{{
!       DO J1=1,NMIN
!          KSUM(J1)=0.0D0
!       ENDDO
!       DO J1=1,NTS
!          IF (PLUS(J1).NE.MINUS(J1)) KSUM(PLUS(J1))=KSUM(PLUS(J1))+EXP(KPLUS(J1)-KMEAN)
!          IF (PLUS(J1).NE.MINUS(J1)) KSUM(MINUS(J1))=KSUM(MINUS(J1))+EXP(KMINUS(J1)-KMEAN)
!       ENDDO
!       DO J1=1,NMIN
!          IF (KSUM(J1).GT.0.0D0) THEN
!             KSUM(J1)=LOG(KSUM(J1))+KMEAN
! !           IF (DEBUG) WRITE(*,'(A,I6,2E20.10)') 'setup> J1,KSUM=',J1,KSUM(J1)
!          ENDIF
!       ENDDO
!       DO J1=1,NTS
! !        IF (DEBUG) WRITE(*,'(A,I6,2E20.10)') 'setup> J1,k+,k-=',J1,KPLUS(J1),KMINUS(J1)
!       ENDDO
!op226 }}}

