      SUBROUTINE CHECKD(X)

      USE COMMONS, ONLY: NATOMS, COMPRESST, PERCOLATET

      IMPLICIT NONE

      INTEGER          :: IVRNO
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), ENERGY, FM, FP, DFA, DFN
      LOGICAL          :: GTEST, STEST, COMPON
      DOUBLE PRECISION, PARAMETER :: ERRLIM = 1.D-06, DELX = 1.D-06
      COMMON /CO/ COMPON

! jwrm2> Turning compression on, if required
      IF (COMPRESST .OR. PERCOLATET) COMPON = .TRUE.

      STEST = .FALSE.

!     Checks gradients

      DO IVRNO = 1, 3*NATOMS

         WRITE(*, *) IVRNO

         GTEST    = .FALSE.
         X(IVRNO) = X(IVRNO) - DELX
         CALL POTENTIAL (X, G, FM, GTEST, STEST)
!         WRITE(*, *) 'Energy minus = ', FM

         X(IVRNO) = X(IVRNO) + 2.D0*DELX
         CALL POTENTIAL (X, G,  FP, GTEST, STEST)
!         WRITE(*, *) 'Energy plus  = ', FP

         GTEST = .TRUE.
         X(IVRNO) = X(IVRNO) - DELX
         CALL POTENTIAL (X, G, ENERGY, GTEST, STEST)
         DFN = (FP - FM) / (2.D0*DELX)
         IF (ABS(DFN) .LT. 1.0D-10) DFN = 0.D0
         DFA = G(IVRNO)

         WRITE(*, *) 'Gradient numerical  = ', DFN
         WRITE(*, *) 'Gradient analytical = ', DFA

         IF (ABS(DFN - DFA) > ERRLIM) WRITE(*, *) IVRNO, DFN, DFA, ABS(DFN-DFA)

      ENDDO

      STOP

      END SUBROUTINE CHECKD
