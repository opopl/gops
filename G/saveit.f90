
      SUBROUTINE GSAVEIT(EREAL,P,NP)

      USE COMMONS
      USE V

      IMPLICIT NONE
      ! sub 
      INTEGER NP
      DOUBLE PRECISION :: EREAL
      DOUBLE PRECISION, DIMENSION(:) ::   P
      ! loc

      INTEGER J1, J2, J3,  NQTOT, NPCALL, CSMIT
      DOUBLE PRECISION AVVAL, CSMRMS

      COMMON /TOT/ NQTOT
      COMMON /PCALL/ NPCALL
!
!  Save the lowest NSAVE distinguishable configurations.
!
!     WRITE(*,'(A,12E15.7)') 'EREAL,ECONV,QMIN',EREAL,ECONV,(QMIN(J1),J1=1,NSAVE)
      DO J1=1,NSAVE
         IF (DABS(EREAL-QMIN(J1)).LT.ECONV) THEN
!
!  These are probably the same - but just to make sure we save the 
!  lowest.
!
            IF (EREAL.LT.QMIN(J1)) THEN
               QMIN(J1)=EREAL
               DO J2=1,3*NATOMS
                  QMINP(J1,J2)=P(J2)
               ENDDO
            ENDIF
            GOTO 10
         ENDIF
         IF (EREAL.LT.QMIN(J1)) THEN

            J2=NSAVE
20          CONTINUE
      
            IF (NSAVE.GT.1) THEN
               QMIN(J2)=QMIN(J2-1)
               FF(J2)=FF(J2-1)
               DO J3=1,3*NATOMS
                  QMINP(J2,J3)=QMINP(J2-1,J3)
               ENDDO

               J2=J2-1
               IF (J2.GE.J1+1) GOTO 20
            ENDIF

            QMIN(J1)=EREAL
            FF(J1)=NQ(NP)
            DO J2=1,3*NATOMS
               QMINP(J1,J2)=P(J2)
            ENDDO

            GOTO 10
         ENDIF
      ENDDO

10    CONTINUE

      DO J1=1,NTARGETS
         IF (EREAL-TARGETS(J1).LT.ECONV) THEN
            IF (NPAR.LT.2) THEN
               WRITE(LFH,'(2(A,I15),A)') 'saveit> Target hit after ',NQTOT,' quenches ',NPCALL,' function calls.'
               WRITE(LFH,'(2(A,F20.10))') 'saveit> Energy=',EREAL,' target=',TARGETS(J1)
            ELSE
               WRITE(LFH,'(A,I1,A,I2,2(A,I15),A)') '[',NP,']saveit> Target hit in parallel run ',NP,&
     &                    ' after ',NQTOT,' quenches ',NPCALL,' function calls.'
               WRITE(LFH,'(A,I1,2(A,F20.10))') '[',NP,']saveit> Energy=',EREAL,' target=',TARGETS(J1)
            ENDIF
            HIT=.TRUE.
         ENDIF
      ENDDO

      RETURN
      END

