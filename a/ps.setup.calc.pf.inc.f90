
!
!  Calculate partition functions for minima. Note that the total partition function
!  is not needed, just the totals for A and B. Since A and B sets are fixed here
!  we don;t need to change the totals.
!
!     PFMEAN=0.0D0
      PFMEAN=-HUGE(1.0D0)
!     PFNORM1=0.0D0 ! use this to calculate ratios without the pe factor
!     PFNORM2=0.0D0 ! use this to calculate ratios with the pe factor

      !op226 {{{

      IF (ENSEMBLE.EQ.'T') THEN
        ! {{{
         IF (TEMPERATURE.LE.0.0D0) THEN
            PRINT '(A)','setup> ERROR - TEMPERATURE=',TEMPERATURE
            STOP
         ENDIF
         DO J1 = 1,NMIN
            PFMIN(J1) = -EMIN(J1)/TEMPERATURE - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))
!           PFMEAN=PFMEAN+PFMIN(J1)
!           PFNORM1=PFNORM1+EXP(- FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1)))
!           PFNORM2=PFNORM2+EXP(PFMIN(J1))
            IF (PFMIN(J1).GT.PFMEAN) PFMEAN=PFMIN(J1)
         ENDDO
         ! }}}
      ELSEIF (ENSEMBLE.EQ.'E') THEN
         ! {{{
         DO J1 = 1,NMIN
            IF (TOTALE.GT.EMIN(J1)) THEN
               PFMIN(J1) = (KAPPA-1)*LOG(TOTALE-EMIN(J1)) - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))
!              PFMEAN=PFMEAN+PFMIN(J1)
!              PFNORM1=PFNORM1+EXP(- FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1)))
!              PFNORM2=PFNORM2+EXP(PFMIN(J1))
               IF (PFMIN(J1).GT.PFMEAN) PFMEAN=PFMIN(J1)
            ELSE
               PFMIN(J1) = -1.0D250
            ENDIF
         ENDDO
         ! }}}
      ELSE
         PRINT*,'ERROR, ENSEMBLE must be set to T or E'
         STOP
      ENDIF

      ! }}}

      IF (DEBUG) THEN
         WRITE(*,'(A,G20.10)') 'setup> mean ln Z=',PFMEAN
!        WRITE(*,'(A)') '     energy        pg order     high T/E prob       Peq'
!        DO J1=1,NMIN
!           WRITE(*,'(F20.10,I6,2G20.10)') EMIN(J1),HORDERMIN(J1), 
!    &                    EXP(-FVIBMIN(J1)/2.0D0-LOG(1.0D0*HORDERMIN(J1))-LOG(PFNORM1)), 
!    &                    EXP(PFMIN(J1)-LOG(PFNORM2))
!        ENDDO
      ENDIF
      DO J1=1,NMIN
         PFMIN(J1) = PFMIN(J1) - PFMEAN
      ENDDO

      PFTOTALB=0.0D0
      DO J1=1,NMINB
         PFTOTALB=PFTOTALB+EXP(PFMIN(LOCATIONB(J1))-PFMIN(LOCATIONB(1)))
      ENDDO
      IF (NMINB.GT.0.0D0) PFTOTALB=LOG(PFTOTALB)+PFMIN(LOCATIONB(1))

      PFTOTALA=0.0D0
      DO J1=1,NMINA
         PFTOTALA=PFTOTALA+EXP(PFMIN(LOCATIONA(J1))-PFMIN(LOCATIONA(1)))
      ENDDO
      IF (NMINA.GT.0.0D0) PFTOTALA=LOG(PFTOTALA)+PFMIN(LOCATIONA(1))



