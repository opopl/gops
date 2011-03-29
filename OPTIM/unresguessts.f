C   OPTIM: A PROGRAM FOR OPTIMIZING GEOMETRIES AND CALCULATING REACTION PATHWAYS
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF OPTIM.
C
C   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
C
C ROUTINE TO GUESS TRANSITION STATES FOR UNRES (DESPITE IT'S NAME...) 
C BY INTERPOLATING BETWEEN DIFFERENT INTERNAL COORDINATE DIHEDRAL VALUES
C DESIGNED TO REPLACE NEB ROUTINE (WHICH IS CALLED FROM CONNECT).
C EQUIVALENT TO DAE'S CHGUESSTS FOR CHARMM!
C
      SUBROUTINE UNRESGUESSTS(Q,ITEST,PTEST,TWISTTYPE,TWISTFRAC,GUESSFAIL,DISTPF)
      USE COMMONS
      USE KEY
      USE MODTWOEND
      USE MODUNRES
      IMPLICIT NONE
C
      DOUBLE PRECISION ANGLE,TWISTFRAC,Q(3*NATOMS)
C
      REAL*8 DIFFPP,SAVEDIFFPP,MAXDIFF2, DUMMYA, RAND, SUMDIFF, DPRAND
C JMC CHANGED DIMENSION OF THE FOLLOWING THREE ARRAYS... WAS MXATMS.
      REAL*8 FINPPSANGLE(4*NRES-9),QPPSANGLE(4*NRES-9),DIFFARRAY(4*NRES-9),DISTPF
      INTEGER IMIN1,IMIN2,IICD,TWISTMODE,TWISTTYPE,NM,NWRONG,I1,J1
      LOGICAL PTEST,ITEST,RANDOM,NORANDOM,GUESSFAIL
      CHARACTER(LEN=18) GUESSFNAME

      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)
      CHARACTER*15 DIHENAME
      DOUBLE PRECISION ENERGY,RMS,GRAD(3*NATOMS)
      INTEGER NWRONGPOL,TWISTMODEPOL,SAVEDIFFPPPOL,NANGLE

      IF (FILTH.EQ.0) THEN
         GUESSFNAME='UNGUESSTS.XYZ'
      ELSE
         WRITE(GUESSFNAME,'(A)') 'UNGUESSTS.XYZ.'//TRIM(ADJUSTL(FILTHSTR))
      ENDIF

C     OPEN(78,FILE='CHGUESSTS.XYZ',STATUS='UNKNOWN')
      OPEN(78,FILE=GUESSFNAME,STATUS='UNKNOWN')
      CALL UNRESDUMP2(Q,78)
C
      DIFFARRAY=0.0D0 ! JMC INITIALISING

      DO I1=1,NRES
         C(1,I1)=FIN(6*(I1-1)+1)
         C(2,I1)=FIN(6*(I1-1)+2)
         C(3,I1)=FIN(6*(I1-1)+3)
         C(1,I1+NRES)=FIN(6*(I1-1)+4)
         C(2,I1+NRES)=FIN(6*(I1-1)+5)
         C(3,I1+NRES)=FIN(6*(I1-1)+6)
C     PRINT *,'FIN IN UNRESGUESSTS: ',FIN(6*(I1-1)+1),FIN(6*(I1-1)+2),FIN(6*(I1-1)+3)
C     PRINT *,'FIN IN UNRESGUESSTS: ',FIN(6*(I1-1)+4),FIN(6*(I1-1)+5),FIN(6*(I1-1)+6)
      END DO
      CALL UPDATEDC
      CALL INT_FROM_CART(.TRUE.,.FALSE.)

      DO I1=1,NRES-3
        FINPPSANGLE(I1)=PHI(I1+3)
      END DO
      DO I1=1,NRES-2
        FINPPSANGLE(I1+NRES-3)=OMEG(I1+1)
C JMC 30/4/03 TRY ADDING BACKBONE AND SIDE CHAIN POLAR ANGLES TO THE INTERPOLATION PROCEDURE...
C THIS SHOULD BE MORE IMPORTANT FOR UNRES THAN FOR CHARMM...
C ORDER IS BB DIHEDRALS, SC DIHEDRALS, BB BOND ANGLES, SC POLARS.
        FINPPSANGLE(I1+2*NRES-5)=THETA(I1+2)
        FINPPSANGLE(I1+3*NRES-7)=ALPH(I1+1)
      END DO

      DO I1=1,NRES
         C(1,I1)=Q(6*(I1-1)+1)
         C(2,I1)=Q(6*(I1-1)+2)
         C(3,I1)=Q(6*(I1-1)+3)
         C(1,I1+NRES)=Q(6*(I1-1)+4)
         C(2,I1+NRES)=Q(6*(I1-1)+5)
         C(3,I1+NRES)=Q(6*(I1-1)+6)
      END DO
      CALL UPDATEDC
      CALL INT_FROM_CART(.TRUE.,.FALSE.)

C USE UNRES GEOMETRY ARRAYS PHI (BB DIHEDRALS) AND OMEG (SC DIHEDRALS)
C NOTE THAT ANGLES ARE IN RADIANS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C TAKE CARE WITH NUMBERING - SEE /UNRES/SRC/READPDB.F (SUBROUTINE INT_FROM_CART)
C FOR SIDE CHAIN DIHEDRALS, THE ACTUAL STORED ARRAYS (ALPHA AND OMEG) CONTAIN ZERO ELEMENTS FOR
C PROPER (I.E. NOT CAPPING) GLYCINES BUT THE VARIABLE ARRAY FROM A CALL TO GEOM_TO_VAR DOES NOT
C CONTAIN THESE ELEMENTS.
C NEED TO REMEMBER NOT TO TRY TO TWIST THEM THOUGH!
C NO ENTRIES IN QPPSANGLE FOR CAPPING 'RESIDUES'.
      DO I1=1,NRES-3
        QPPSANGLE(I1)=PHI(I1+3)
      END DO
      DO I1=1,NRES-2
        QPPSANGLE(I1+NRES-3)=OMEG(I1+1)
        QPPSANGLE(I1+2*NRES-5)=THETA(I1+2)
        QPPSANGLE(I1+3*NRES-7)=ALPH(I1+1)
      END DO
C JMC NOTE THAT THE Q INTERNAL COORD SET IS NOW SAVED IN THE UNRES INT COOR COMMON BLOCK...
C PUT THE TS GUESS COORDS INTO COMMON BLOCK BEFORE EXITING THIS SUBROUTINE.

C NOW DECIDE WHICH PHI/PSI OR SIDECHAIN ANGLE TO TWIST
C
C BASED ON TWISTTYPE
C TWISTTYPE = 1  MEANS TAKE ONE WITH BIGGEST DIFFERENCE AND INTERPOLATE BETWEEN
C      THE TWO VALUES USING TWISTFRAC AS THE FRACTION
C TWISTTYPE = 2; INTERPOLATES LIKE 1 BUT SETS CHOSEN ANGLE TO THE NEAREST OF -120, 0, 120 DEGREES
C I.E. MAXIMA OF THE DIHEDRAL POTENTIAL (WHICH IS K(1+COS(3*PHI)) FOR PHI AND PSI ANGLES. ! CHARMM
C IN FACT K = 0 FOR PSI, SO THIS METHOD MAY BE A BIT FUTILE FOR PSI ANGLES, BUT IT MAY GIVE SENSIBLE ! CHARMM
C GEOMETRIES ANYWAY) ! CHARMM
C
C TWISTTYPE =3; LIKE 1, BUT ALSO INTERPOLATES THE DIHEDRAL EITHER SIDE OF THE MAXIMUM
C
C TWISTTYPE =4; DOES ON ONE DIHEDRAL, CHOSEN WITH PROBABILITY BASED ON SIZE OF DISPLACEMENT
C
C TWISTTYPE =5; IF ONLY ONE DIHEDRAL DIFFERS BY >60DEG THEN INTERPOLATES ON ONE DIHEDRAL,
C               IF MORE THAN ONES DIFFERS THEN PROCEEDS LIKE RANDOM MODE (TWISTTYPE =4)
C
C TWISTTYPE =6; IF ONLY ONE DIHEDRAL DIFFERS BY >60DEG THEN INTERPOLATES THAT DIHEDRAL,
C               AND THE DIHEDRALS EITHER SIDE (LIKE TT=3)
C               IF MORE THAN ONES DIFFERS THEN PROCEEDS LIKE RANDOM MODE (TWISTTYPE =4)
C
C TWISTTYPE =7; JUST INTERPOLATE ALL DIHEDRALS
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C JMC NOTE FOR NOW, ONLY THESE OPTIONS BELOW ARE WORKING!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C TWISTTYPE =8; INTERPOLATE ALL BACKBONE ANGLES
C TWISTTYPE =9; INTERPOLATE LARGEST DIHEDRAL AND LARGEST POLAR ANGLE
C TWISTTYPE =10; JUST INTERPOLATE ALL ANGLES
C
      RANDOM=.FALSE.
      MAXDIFF2 = 0.0D0
      NWRONG=0
      NWRONGPOL=0
C
C TURN RANDOM DISPLACEMENTS OFF ONCE TWO MINIMA ARE CLOSE ENOUGH TOGETHER FOR NEB
C TO BE SUCCESSFUL
C
C     IF (DISTPF.LT.RANDOMCUTOFF) THEN
C        NORANDOM=.TRUE.
C     ELSE
         NORANDOM=.FALSE.
C     ENDIF

      DO I1=1,NPHI+NRES-2
         IF (I1.GT.NPHI) THEN
            IF (ITYPE(I1-NPHI+1).EQ.10) GOTO 100 ! GLYCINE
         END IF
C        WRITE(*,'(A,I6,2F15.10)') 'FINS QS',I1,FINPPSANGLE(I1),QPPSANGLE(I1)
         DIFFPP = FINPPSANGLE(I1) - QPPSANGLE(I1)
C
C NEXT TWO LINES ARE MEANT TO ENSURE THAT YOU ALWAYS INTERPOLATE
C ALONG THE SHORTEST DISTANCE BETWEEN THE DIHEDRAL ANGLES.
C
         IF (DIFFPP.GT.PI) DIFFPP = DIFFPP-2.0D0*PI
         IF (DIFFPP.LT.-PI) DIFFPP = DIFFPP+2.0D0*PI

         DIFFARRAY(I1)=DIFFPP

         IF ((DIFFPP*DIFFPP).GT.MAXDIFF2) THEN
            MAXDIFF2=DIFFPP*DIFFPP
            TWISTMODE=I1
            SAVEDIFFPP=DIFFPP
         ENDIF

C JMC         IF (ABS(DIFFPP).GT.60.0D0) NWRONG=NWRONG+1
         IF (ABS(DIFFPP).GT.PI/3.0D0) NWRONG=NWRONG+1

100   CONTINUE 
      ENDDO

C JMC DON'T DUPLICATE WORK FROM ABOVE DO LOOP...
C REMEMBER THE POLAR ANGLES RUN FROM 0 TO PI, WHEREAS DIHEDRALS GO FROM -PI TO PI.
      DO I1=NPHI+NRES-1,NPHI+NTHETA+2*NRES-4
         IF (I1.GT.NPHI+NRES-2+NTHETA) THEN
            IF (ITYPE(I1-NPHI-NRES+2-NTHETA+1).EQ.10) GOTO 200 ! GLYCINE
         END IF
C        WRITE(*,'(A,I6,2F15.10)') 'FINS QS',I1,FINPPSANGLE(I1),QPPSANGLE(I1)
         DIFFPP = FINPPSANGLE(I1) - QPPSANGLE(I1)
         DIFFARRAY(I1)=DIFFPP
         IF ((DIFFPP*DIFFPP).GT.MAXDIFF2) THEN
            MAXDIFF2=DIFFPP*DIFFPP
            TWISTMODEPOL=I1
            SAVEDIFFPPPOL=DIFFPP
         ENDIF

C JMC         IF (ABS(DIFFPP).GT.60.0D0) NWRONG=NWRONG+1
         IF (ABS(DIFFPP).GT.PI/3.0D0) NWRONGPOL=NWRONGPOL+1

200   CONTINUE
      ENDDO

C
C NOW DO TWISTING 
C

      IF (TWISTTYPE.EQ.7) THEN
         DO I1=1,NPHI
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
C JMC            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-DIFFARRAY(I1))
            IF (TWISTFRAC.LT.0.D0) THEN
               ANGLE = 0.5D0*(2.0D0*PI-ABS(DIFFARRAY(I1))) ! SO ANGLE WILL ALWAYS BE BETWEEN PI/2 AND PI (ALWAYS > 0)
               IF (DIFFARRAY(I1).GT.0.0D0) ANGLE = -ANGLE ! JMC NEED TO TEST THIS!! OR DO WE NEED -1.0D0*ANGLE??
            ENDIF
            PHI(I1+3)=PHI(I1+3)+ANGLE
         ENDDO
         DO I1=NPHI+1,NPHI+NRES-2
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
            IF (TWISTFRAC.LT.0.D0) THEN
               ANGLE = 0.5D0*(2.0D0*PI-ABS(DIFFARRAY(I1))) ! SO ANGLE WILL ALWAYS BE BETWEEN PI/2 AND PI (ALWAYS > 0)
               IF (DIFFARRAY(I1).GT.0.0D0) ANGLE = -ANGLE ! JMC NEED TO TEST THIS!! OR DO WE NEED -1.0D0*ANGLE??
            ENDIF
C JMC            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-DIFFARRAY(I1))
            OMEG(I1+1-NPHI)=OMEG(I1+1-NPHI)+ANGLE
         ENDDO
         GOTO 20
      ENDIF

      IF (TWISTTYPE.EQ.10) THEN
         DO I1=1,NPHI
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
C JMC            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-DIFFARRAY(I1))
            IF (TWISTFRAC.LT.0.D0) THEN
               ANGLE = 0.5D0*(2.0D0*PI-ABS(DIFFARRAY(I1))) ! SO ANGLE WILL ALWAYS BE BETWEEN PI/2 AND PI (ALWAYS > 0)
               IF (DIFFARRAY(I1).GT.0.0D0) ANGLE = -ANGLE ! JMC NEED TO TEST THIS!! OR DO WE NEED -1.0D0*ANGLE??
            ENDIF
            PHI(I1+3)=PHI(I1+3)+ANGLE
         ENDDO
         DO I1=NPHI+1,NPHI+NRES-2
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
C JMC            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-DIFFARRAY(I1))
            IF (TWISTFRAC.LT.0.D0) THEN
               ANGLE = 0.5D0*(2.0D0*PI-ABS(DIFFARRAY(I1))) ! SO ANGLE WILL ALWAYS BE BETWEEN PI/2 AND PI (ALWAYS > 0)
               IF (DIFFARRAY(I1).GT.0.0D0) ANGLE = -ANGLE ! JMC NEED TO TEST THIS!! OR DO WE NEED -1.0D0*ANGLE??
            ENDIF
            OMEG(I1+1-NPHI)=OMEG(I1+1-NPHI)+ANGLE
         ENDDO
         DO I1=NPHI+NRES-1,NPHI+NRES-2+NTHETA
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*DIFFARRAY(I1) ! 'GOING THE LONG WAY ROUND' DOESN'T APPLY FOR BOND ANGLES
            THETA(I1-NPHI-NRES+4)=THETA(I1-NPHI-NRES+4)+ANGLE
         ENDDO
         DO I1=NPHI+NRES-1+NTHETA,NPHI+2*NRES-4+NTHETA
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*DIFFARRAY(I1) ! 'GOING THE LONG WAY ROUND' DOESN'T APPLY FOR BOND ANGLES
            ALPH(I1-NPHI-NTHETA-NRES+3)=ALPH(I1-NPHI-NTHETA-NRES+3)+ANGLE
         ENDDO
         GOTO 20
      ENDIF

      IF (TWISTTYPE.EQ.8) THEN
C JMC BACKBONE ANGLES ONLY
         DO I1=1,NPHI
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
C JMC            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-DIFFARRAY(I1))
            IF (TWISTFRAC.LT.0.D0) THEN
               ANGLE = 0.5D0*(2.0D0*PI-ABS(DIFFARRAY(I1))) ! SO ANGLE WILL ALWAYS BE BETWEEN PI/2 AND PI (ALWAYS > 0)
               IF (DIFFARRAY(I1).GT.0.0D0) ANGLE = -ANGLE ! JMC NEED TO TEST THIS!! OR DO WE NEED -1.0D0*ANGLE??
            ENDIF
            PHI(I1+3)=PHI(I1+3)+ANGLE
         ENDDO
         DO I1=NPHI+NRES-1,NPHI+NRES-2+NTHETA
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*DIFFARRAY(I1) ! 'GOING THE LONG WAY ROUND' DOESN'T APPLY FOR BOND ANGLES
            THETA(I1-NPHI-NRES+4)=THETA(I1-NPHI-NRES+4)+ANGLE
         ENDDO
         GOTO 20
      ENDIF

      IF (TWISTTYPE.EQ.9) THEN
            ANGLE=TWISTFRAC*DIFFARRAY(TWISTMODE)
C JMC            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-DIFFARRAY(TWISTMODE))
            IF (TWISTFRAC.LT.0.D0) THEN
               ANGLE = 0.5D0*(2.0D0*PI-ABS(DIFFARRAY(TWISTMODE))) ! SO ANGLE WILL ALWAYS BE BETWEEN PI/2 AND PI (ALWAYS > 0)
               IF (DIFFARRAY(TWISTMODE).GT.0.0D0) ANGLE = -ANGLE ! JMC NEED TO TEST THIS!! OR DO WE NEED -1.0D0*ANGLE??
            ENDIF
            PHI(TWISTMODE+3)=PHI(TWISTMODE+3)+ANGLE

            ANGLE=TWISTFRAC*DIFFARRAY(TWISTMODE+NPHI+NRES-2)
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*DIFFARRAY(TWISTMODE+NPHI+NRES-2) ! 'GOING THE LONG WAY ROUND' DOESN'T APPLY FOR BOND ANGLES
            THETA(TWISTMODE-NPHI-NRES+4)=THETA(TWISTMODE-NPHI-NRES+4)+ANGLE
            ANGLE=TWISTFRAC*DIFFARRAY(TWISTMODE+1+NPHI+NRES-2)
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*DIFFARRAY(TWISTMODE+1+NPHI+NRES-2) ! 'GOING THE LONG WAY ROUND' DOESN'T APPLY FOR BOND ANGLES
            THETA(TWISTMODE+1-NPHI-NRES+4)=THETA(TWISTMODE+1-NPHI-NRES+4)+ANGLE
         GOTO 20
      ENDIF

      IF ((TWISTTYPE.EQ.5).OR.(TWISTTYPE.EQ.6)) THEN
         IF (NWRONG.GT.2) THEN
C            WRITE (*,'(A)') 'MORE THAN ONE DIHEDRAL DISPLACED - UNLIKELY TO BE A DIRECT CONNECTION'
            WRITE (*,'(A)') 'MORE THAN TWO DIHEDRALS DISPLACED - UNLIKELY TO BE A DIRECT CONNECTION'
            IF (NORANDOM) THEN
               WRITE (*,'(A)') 'SWITCHING TO NEB'
               GUESSFAIL=.TRUE.
               RETURN
            ELSE
               WRITE (*,'(A)') 'CHOOSING A MODE TO TWIST AT RANDOM'
               RANDOM=.TRUE.
C              STOP
            ENDIF
         ENDIF
      ENDIF

      ANGLE=TWISTFRAC*SAVEDIFFPP
C JMC WHAT IF SAVEDIFFPP IS LT 0?
      IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-SAVEDIFFPP) ! THIS WILL NOT WORK PROPERLY - DON'T USE!!!
      IF (TWISTTYPE.EQ.2) THEN
         DUMMYA=QPPSANGLE(TWISTMODE)+ANGLE
C JMC         IF ((DUMMYA.GT.-180.0D0).AND.(DUMMYA.LT.-60.0D0)) DUMMYA=-120.0D0
C JMC         IF ((DUMMYA.GT.-60.0D0).AND.(DUMMYA.LT.60.0D0)) DUMMYA=0.0D0
C JMC         IF ((DUMMYA.GT.60.0D0).AND.(DUMMYA.LT.180.0D0)) DUMMYA=120.0D0
         IF ((DUMMYA.GT.-PI).AND.(DUMMYA.LT.-PI/3.0D0)) DUMMYA=-2.0D0*PI/3.0D0
         IF ((DUMMYA.GE.-PI/3.0D0).AND.(DUMMYA.LT.PI/3.0D0)) DUMMYA=0.0D0
         IF ((DUMMYA.GE.PI/3.0D0).AND.(DUMMYA.LE.PI)) DUMMYA=2.0D0*PI/3.0D0
         ANGLE=DUMMYA-QPPSANGLE(TWISTMODE)
      ENDIF
      
      IF ((TWISTTYPE.EQ.4).OR.RANDOM) THEN
         SUMDIFF=0.D0
         DO I1=1,NPHI+NRES-2
            SUMDIFF=SUMDIFF+ABS(DIFFARRAY(I1))
         ENDDO
         RAND=DPRAND()*SUMDIFF
         PRINT *,'RAND',RAND
         SUMDIFF=0.D0
         DO I1=1,NPHI+NRES-2
            SUMDIFF=SUMDIFF+ABS(DIFFARRAY(I1))
C              PRINT *,'DIFFARRAY ',DIFFARRAY(I1)
            IF (SUMDIFF.GT.RAND) THEN 
               TWISTMODE=I1
               ANGLE=TWISTFRAC*DIFFARRAY(I1)
               PRINT *,'ANGLE ',ANGLE
               IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5*(2.0D0*PI-SAVEDIFFPP)
               GOTO 10
            ENDIF
         ENDDO
10      CONTINUE
C JMC HUH?
C       IF (RANDOM) THEN
C JMC           ANGLE=DPRAND()*60.D0
C          ANGLE=DPRAND()*PI/3.0D0
C          IF (DIFFARRAY(TWISTMODE).LT.0.D0) ANGLE=-1.D0*ANGLE
C       ENDIF
      ENDIF

C JMC      WRITE (*,'(A20,I3,A2,1X,F10.5,1X,A8)') 'TWISTING PHI/PSI DIHEDRAL ',TWISTMODE,' BY ',ANGLE,' DEGREES'
      WRITE (*,'(A20,I3,A2,1X,F10.5,1X,A8)') 'TWISTING PHI/PSI DIHEDRAL ',TWISTMODE,' BY ',ANGLE,' RADIANS'

      IF (TWISTMODE.LE.NPHI) THEN
         PHI(TWISTMODE+3)=PHI(TWISTMODE+3)+ANGLE
      ELSE
         OMEG(TWISTMODE+1-NPHI)=OMEG(TWISTMODE+1-NPHI)+ANGLE
      END IF

      IF ((TWISTTYPE.EQ.3).OR.((TWISTTYPE.EQ.6).AND.(.NOT.RANDOM))) THEN
         NM=TWISTMODE-1
         IF (NM.GE.1) THEN
            DIFFPP = FINPPSANGLE(NM) - QPPSANGLE(NM)
C JMC            IF (DIFFPP.GT.180.0) DIFFPP = DIFFPP-360.D0
C JMC            IF (DIFFPP.GT.180.0) DIFFPP = DIFFPP-360.D0
            IF (DIFFPP.LT.-PI) DIFFPP = DIFFPP+2.0D0*PI
            IF (DIFFPP.GT.PI) DIFFPP = DIFFPP-2.0D0*PI
            ANGLE=TWISTFRAC*DIFFPP
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-SAVEDIFFPP)
C JMC            IICD=PHIPSI(NM)
C JMC            CALL TWISTCH(IICD,ANGLE)
            PHI(NM+3)=PHI(NM+3)+ANGLE
         ENDIF
         NM=TWISTMODE+1
C JMC         IF (NM.LE.NPHIPSI) THEN
         IF (NM.LE.NPHI) THEN
            DIFFPP = FINPPSANGLE(NM) - QPPSANGLE(NM)
            IF (DIFFPP.GT.PI) DIFFPP = DIFFPP-2.0D0*PI
            IF (DIFFPP.LT.-PI) DIFFPP = DIFFPP+2.0D0*PI
            ANGLE=TWISTFRAC*DIFFPP
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-SAVEDIFFPP)
C JMC            IICD=PHIPSI(NM)
C JMC            CALL TWISTCH(IICD,ANGLE)
            PHI(NM+3)=PHI(NM+3)+ANGLE
         ENDIF
       ENDIF
C
20    CONTINUE

      CALL CHAINBUILD

      DO J1=1,NRES
         Q(6*(J1-1)+1)=C(1,J1)
         Q(6*(J1-1)+2)=C(2,J1)
         Q(6*(J1-1)+3)=C(3,J1)
         Q(6*(J1-1)+4)=C(1,J1+NRES)
         Q(6*(J1-1)+5)=C(2,J1+NRES)
         Q(6*(J1-1)+6)=C(3,J1+NRES)
      END DO

      CALL UNRESDUMP2(Q,78)
      CALL UNRESDUMP2(FIN,78)

      CLOSE(78)

C JMC TESTING
C     CALL POTENTIAL(Q,ENERGY,GRAD,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
C     PRINT *,'ENERGY,RMS IN CHGUESSTS ',ENERGY,RMS
C     DIHENAME='QTS'
C     CALL PRINTDIHE(Q,QPPSANGLE,NANGLE,DIHENAME)

      RETURN

      END

      SUBROUTINE UNRESGUESSMIN(Q,PTEST,TWISTTYPE,NGUESS)

      USE COMMONS
      USE MODTWOEND
      USE MODUNRES
      USE KEY
      USE VARS
      IMPLICIT NONE

      DOUBLE PRECISION Q(3*NATOMS),TWISTFRAC,DISTPF,TMPQ(3*NATOMS)
      LOGICAL PTEST,GUESSFAIL
      INTEGER TWISTTYPE,I1,NGUESS,J1,K1
      DOUBLE PRECISION Q1(3*NATOMS),Q2(3*NATOMS),Q3(3*NATOMS)
C     COMMON /MINARRAY/ MYQMINSAVE(3*NATOMS,100),COUNTER   ! NOTE ARBITRARILY CHOOSING 100, COULD HAVE AS A PARAMETER BUT CAN'T BE BOTHERED NOW...
      CHARACTER*5 ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE

      CALL NEWMINDIST(Q,FIN,NATOMS,DISTPF,.FALSE.,.FALSE.,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
      PRINT *,'Q, FIN ',DISTPF

      TMPQ=Q

      COUNTER=0
      DO I1=1,NGUESS
         TWISTFRAC=1.0D0*I1/(NGUESS+1)
         CALL UNRESGUESSTS(Q,.FALSE.,PTEST,TWISTTYPE,TWISTFRAC,GUESSFAIL,DISTPF)
         Q=TMPQ
         PRINT *,'GUESSFAIL ',GUESSFAIL
      END DO

C NOW DO SOME FUNKY MIND STUFF ON MYQMINSAVE...
      DO I1=1,COUNTER
         DO K1=1,3*NATOMS
            Q1(K1)=MYQMINSAVE(K1,I1)
            Q2(K1)=TMPQ(K1)
            Q3(K1)=FIN(K1)
         END DO
         CALL NEWMINDIST(Q1,Q2,NATOMS,DISTPF,.FALSE.,.FALSE.,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
         PRINT *,'DISTPF ',DISTPF,I1,' START'
         CALL NEWMINDIST(Q1,Q3,NATOMS,DISTPF,.FALSE.,.FALSE.,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
         PRINT *,'DISTPF ',DISTPF,I1,' FIN'
         DO J1=I1+1,COUNTER
            DO K1=1,3*NATOMS
               Q1(K1)=MYQMINSAVE(K1,I1)
               Q2(K1)=MYQMINSAVE(K1,J1)
            END DO
            CALL NEWMINDIST(Q1,Q2,NATOMS,DISTPF,.FALSE.,.FALSE.,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
            PRINT *,'DISTPF ',DISTPF,I1,J1
         END DO
      END DO

      STOP

      END SUBROUTINE UNRESGUESSMIN
