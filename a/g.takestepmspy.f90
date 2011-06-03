      SUBROUTINE TAKESTEPMSPY (NP)

!     THIS ROUTINE TAKES STEP FOR MULTI-SITE RIGID BODIES ENSURING NO OVERLAP

      USE COMMONS  

      IMPLICIT NONE

      INTEGER          :: NP, JMAX, JMAX2, REALNATOMS, OFFSET, PTINDX
      INTEGER          :: J1, J2, J3, J4, J5, J6, I, J
      LOGICAL          :: OVRLPT

      DOUBLE PRECISION :: PI, DUMMY, DUMMY2
      DOUBLE PRECISION :: DIST(3*NATOMS/2), XMASS, YMASS, ZMASS, DMAX, VMAX, VMAX2
      DOUBLE PRECISION :: VMIN, CMMAX, CMDIST(NATOMS/2), LOCALSTEP
      DOUBLE PRECISION :: DPRAND, RANDOM, THETA, PHI, THETA2, CT, ST, XMIN, FMIN
      DOUBLE PRECISION :: EZRI1(3,3), EZRI2(3,3), EZRJ1(3,3), EZRJ2(3,3)
      DOUBLE PRECISION :: PST(NPYSITE,3), RMI(3,3), RMJ(3,3), E(3,3), AE(3,3), BE(3,3), I3(3,3) 
      DOUBLE PRECISION :: RI(3), RJ(3), RISITE(3), RJSITE(3), RIJ(3)
      DOUBLE PRECISION :: P(3), ABSRIJ, RCUT, ECFVAL

      PI         = 4.D0*ATAN(1.0D0)
      RCUT       = 2.D0 * MAX(MAXVAL(PYA1),MAXVAL(PYA2))
      REALNATOMS = NATOMS/2
      OFFSET     = 3 * REALNATOMS
      I3(:,:)    = 0.D0

      DO I = 1, 3

         I3(I,I) = 1.D0

      ENDDO

      CALL DEFMSPY(PST)

      DO J1 = 1,REALNATOMS

         J2     = 3*J1
         DUMMY2 = COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2
         IF (DUMMY2 .GT. RADIUS) THEN
            WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,R**2,x,y,z:', J1, RADIUS, DUMMY2, COORDS(J2-2,NP), &
                                       COORDS(J2-1,NP), COORDS(J2,NP)
            PRINT*, 'initial coordinate outside container -- increase container radius'
            STOP
         END IF

      END DO

      DO J1 = 1,3*NATOMS
         COORDSO(J1,NP) = COORDS(J1,NP)
      END DO

      DO J1 = 1,NATOMS/2
         VATO(J1,NP) = VAT(J1,NP)
      END DO

!     FIND THE CENTRE OF MASS

      XMASS = 0.0D0
      YMASS = 0.0D0
      ZMASS = 0.0D0

      DO J1 = 1,NATOMS/2

         XMASS = XMASS + COORDS(3*(J1-1)+1,NP)
         YMASS = YMASS + COORDS(3*(J1-1)+2,NP)
         ZMASS = ZMASS + COORDS(3*(J1-1)+3,NP)

      ENDDO

      XMASS = XMASS/(REALNATOMS)
      YMASS = YMASS/(REALNATOMS)
      ZMASS = ZMASS/(REALNATOMS)

!     Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
!     and the pair energy of the most tightly bound atom, VMIN. An angular step is
!     taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
!     DMAX (or CMMAX from CM of the cluster).

      DMAX  =  1.0D0
      VMAX  = -1.0D3
      VMAX2 = -1.0D3
      VMIN  =  0.0D0
      CMMAX =  1.0D0

      DO J1 = 1, REALNATOMS

         J2 = 3*J1
         DIST(J1)   = DSQRT( COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2)
         CMDIST(J1) = DSQRT((COORDS(J2-2,NP)-XMASS)**2+(COORDS(J2-1,NP)-YMASS)**2+(COORDS(J2,NP)-ZMASS)**2)
         IF (CMDIST(J1) .GT. CMMAX) CMMAX = CMDIST(J1)
         IF (DIST(J1) .GT. DMAX) DMAX = DIST(J1)
         IF (VAT(J1,NP) .GT. VMAX) THEN
            VMAX = VAT(J1,NP)
            JMAX = J1
         ELSE IF ((VAT(J1,NP).LT. VMAX) .AND. (VAT(J1,NP) .GT. VMAX2)) THEN
              VMAX2 = VAT(J1,NP)
              JMAX2 = J1
         ENDIF
         IF (VAT(J1,NP) .LT. VMIN) VMIN = VAT(J1,NP)

      ENDDO

      IF (VAT(JMAX,NP) > (ASTEP(NP)*VMIN) .AND. (.NOT.NORESET)) THEN

         J2 = 3*JMAX
         THETA           = DPRAND()*PI
         PHI             = DPRAND()*PI*2.0D0
         COORDS(J2-2,NP) = XMASS + (CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
         COORDS(J2-1,NP) = YMASS + (CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
         COORDS(J2,NP)   = ZMASS + (CMMAX+1.0D0)*DCOS(THETA)
         DUMMY           = COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2

         IF (DUMMY > RADIUS) THEN

            DUMMY           = DSQRT(RADIUS*0.99D0/DUMMY)
            COORDS(J2-2,NP) = COORDS(J2-2,NP)*DUMMY
            COORDS(J2-1,NP) = COORDS(J2-1,NP)*DUMMY
            COORDS(J2,NP)   = COORDS(J2,NP)*DUMMY

         END IF

      ENDIF

      DO J1 = 1, NATOMS

         J3 = 3*J1

         LOCALSTEP = STEP(NP)

         IF (J1 > REALNATOMS) THEN

            LOCALSTEP = 0.0D0
            IF (OMOVE(NP)) LOCALSTEP = OSTEP(NP)

         ELSE IF (J1 <= REALNATOMS) THEN

            LOCALSTEP = 0.0D0
            IF (TMOVE(NP)) LOCALSTEP = STEP(NP)

         END IF

!     CHECK FOR OVERLAP

         OVRLPT = .TRUE.

         DO WHILE (OVRLPT)

            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3-2,NP) = COORDS(J3-2,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3-1,NP) = COORDS(J3-1,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3,NP)   = COORDS(J3,NP) + LOCALSTEP*RANDOM
          
            IF (J1 <= REALNATOMS) THEN

               PTINDX = J1
               J5     = OFFSET + J3
               RI(:)  = COORDS(J3-2:J3,NP)
               P(:)   = COORDS(J5-2:J5,NP)

            ELSE 

               PTINDX = J1 - REALNATOMS
               J5     = J3 - OFFSET
               RI(:)  = COORDS(J5-2:J5,NP)
               P(:)   = COORDS(J3-2:J3,NP)

            ENDIF  

!     ROTATION MATRIX

            THETA  = DSQRT(DOT_PRODUCT(P,P))

            IF (THETA == 0.D0) THEN

               RMI = I3

            ELSE

               THETA2 = THETA * THETA
               CT      = COS(THETA)
               ST      = SIN(THETA)
               E(:,:)  = 0.D0
               E(1,2)  = -P(3)
               E(1,3)  =  P(2)
               E(2,3)  = -P(1)
               E(2,1)  = -E(1,2)
               E(3,1)  = -E(1,3)
               E(3,2)  = -E(2,3)
               E       = E/THETA

               RMI     = I3 + (1.D0-CT)*MATMUL(E,E) + E*ST
           
            ENDIF

            DO I = 1, NPYSITE

!     OBTAIN THE SITE POSITION IN THE SPACE-FIXED FRAME  

               RISITE = RI +  MATMUL(RMI,PST(I,:))
                  
               DO J2 = 1, REALNATOMS

                  IF (J2 == PTINDX) CYCLE

                  J4    = 3*J2
                  J6    = OFFSET + J4
                  RJ(:) = COORDS(J4-2:J4,NP)
                  P(:)  = COORDS(J6-2:J6,NP)

!     ROTATION MATRIX

                  THETA  = DSQRT(DOT_PRODUCT(P,P))

                  IF (THETA == 0.D0) THEN

                     RMJ = I3

                  ELSE
 
                     THETA2 = THETA * THETA
                     CT      = COS(THETA)
                     ST      = SIN(THETA)
                     E(:,:)  = 0.D0
                     E(1,2)  = -P(3)
                     E(1,3)  =  P(2)
                     E(2,3)  = -P(1)
                     E(2,1)  = -E(1,2)
                     E(3,1)  = -E(1,3)
                     E(3,2)  = -E(2,3)
                     E       = E/THETA

                     RMJ     = I3 + (1.D0-CT)*MATMUL(E,E) + E*ST

                  ENDIF

                  DO J = 1, NPYSITE

!     OBTAIN THE SITE POSITION IN THE SPACE-FIXED FRAME   

                     RJSITE = RJ + MATMUL(RMJ,PST(J,:))

                     RIJ    = RISITE - RJSITE
                     ABSRIJ = DSQRT(DOT_PRODUCT(RIJ,RIJ))
                       
                     IF (ABSRIJ < RCUT) THEN 

!     DETERMINE ELLIPTIC CONTACT FUNCTION

                        CALL SITEBF(I, EZRI1, EZRI2) 

                        AE = MATMUL(RMI,(MATMUL(EZRI1(:,:),(TRANSPOSE(RMI)))))
                  
                        CALL SITEBF(J, EZRJ1, EZRJ2)

                        BE = MATMUL(RMJ,(MATMUL(EZRJ1(:,:),(TRANSPOSE(RMJ)))))

                        CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE, BE, RIJ, XMIN, FMIN)

                        ECFVAL = - FMIN
  
                        IF (ECFVAL >= 1.D0) THEN
                           OVRLPT = .FALSE.
                        ENDIF

                     ELSE

                        OVRLPT = .FALSE.

                     ENDIF

                  ENDDO  ! END LOOP OVER J

               ENDDO  ! END LOOP OVER I

            ENDDO  ! END LOOP OVER J2

         ENDDO  ! END WHILE
 
      ENDDO  ! END LOOP OVER J1   

      END SUBROUTINE TAKESTEPMSPY
