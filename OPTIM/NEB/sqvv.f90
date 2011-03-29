!   NEB MODULE IS AN IMPLEMENTATION OF THE NUDGED ELASTIC BAND METHOD FOR PERFORMING DOUBLE-ENDED PATHWAY SEARCHES.
!   COPYRIGHT (C) 2003-2006 SEMEN A. TRYGUBENKO AND DAVID J. WALES
!   THIS FILE IS PART OF NEB MODULE. NEB MODULE IS PART OF OPTIM.
!
!   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
!   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
!   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
!   (AT YOUR OPTION) ANY LATER VERSION.
!
!   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
!   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
!   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
!   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
!
!   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
!   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
!   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
!
MODULE MINIMISER2
     IMPLICIT NONE
     CONTAINS

     ! THIS WORK WAS DESCRIBED IN DETAIL IN: S. A. TRYGUBENKO AND D. J. WALES, `A DOUBLY NUDGED ELASTIC BAND METHOD FOR FINDING
     ! TRANSITION STATES', J. CHEM. PHYS., 120, 2082-2094 (2004). SUMMARY IS AVAILABLE ONLINE AT
     ! HTTP://WWW-WALES.CH.CAM.AC.UK/~SAT39/DNEBTESTS/
     SUBROUTINE NEBSQVV(N)
          USE KEYMINIMIZER
          USE NEBDATA
          USE KEYSQVV
          USE KEYNEB,ONLY: NIMAGE,MOREPRINTING,DUMPNEBXYZ,DUMPNEBXYZFREQ,DUMPNEBEOS,DUMPNEBEOSFREQ,DEBUG,&
          &DUMPNEBPTS,DUMPNEBPTSFREQ
          USE GRADIENTS
          USE NEBUTILS
          USE CHARUTILS
          IMPLICIT NONE

          INTEGER,INTENT(IN) :: N

          INTEGER(4) :: I,Q,IMAX
          DOUBLE PRECISION :: PROJ,GNORM
          DOUBLE PRECISION,PARAMETER   :: MASS=1.0D0
          DOUBLE PRECISION,ALLOCATABLE :: V(:),DX(:)
          CHARACTER(LEN=1) :: QUENCHING,QVVTYPE

          ! WHERE TO QUENCH VELOCITY?
          QVVTYPE = "2" ! 2 == SQVV

          ALLOCATE(V(N),DX(N))
          IF (MOREPRINTING) THEN
               PRINT '(/1X,A)', REPEAT(">",21)//" ENTERING SQVV MINIMISER "//REPEAT(">",21)
               WRITE(*,'(1X,A)') REPEAT("_",91)
               WRITE(*,'(A6,A20,A20,A9,2A8,2A9)') 'ITER','ENERGY','RMS FORCE','AV.DEV','PATH','STEP','VEL'
               WRITE(*,'(1X,A)') REPEAT("_",91)
          ENDIF

          IF (DEBUG) CALL DUMPFILES("B")
          
          IF (SQVVGUESS) THEN
               IMAX = NITERSQVVGUESSMAX
          ELSE
               IMAX = NITERMAX
          ENDIF

          DO NITERDONE=1,IMAX
               CALL NEBGRADIENT
               IF (QVVTYPE == "4") CALL QUENCH
               IF (BADTAU) EXIT
               
               ! SECOND VELOCITY UPDATE
               IF (READVEL) THEN
                    PRINT *, "VELOCITY READ FROM FILE VEL."
                    ! READ VELOCITY VECTOR FROM PREVIOUS ITERATION, ASSUME
                    ! GEOMETRY WAS READ AS WELL ALREADY
                    OPEN(UNIT=333,FILE="VEL",STATUS="OLD",FORM="UNFORMATTED")
                    READ(UNIT=333) V
                    CLOSE(333)
                    V = V - (0.5D0*DT/MASS)*G
               ELSE
                    IF (NITERDONE==1) THEN
                         V = 0.0D0
                    ELSE
                         V = V - (0.5D0*DT/MASS)*G
                    ENDIF
               ENDIF

                    IF (QVVTYPE == "1") CALL QUENCH
                    
               ! CALCULATE THE STEP AND LIMIT IT IF APPROPRIATE
               DX = V*DT - (0.5D0*DT**2/MASS)*G ! ZSUNYLYSYA PO KOORDYNATAH; V CHASI - NA DT
               IF (LIMITSQVVSTEP) WHERE (ABS(DX) > STEPDOFMAX) DX = STEPDOFMAX*DX/ABS(DX)
               
               ! CALCULATE THE STEP PER IMAGE USING DX
               DO I=1,NIMAGE
                    STEPIMAGE(I) = SQRT(DOT_PRODUCT(DX(NOPT*(I-1)+1:NOPT*I),DX(NOPT*(I-1)+1:NOPT*I)))
               ENDDO
               STEPTOT = SUM(STEPIMAGE)/NIMAGE
               
               X = X + DX ! COORDINATE UPDATE
                    IF (QVVTYPE == "2") CALL QUENCH
               V =     V    - (0.5D0*DT  / MASS)*G ! FIRST VELOCITY UPDATE
                    IF (QVVTYPE == "3") CALL QUENCH
               
               CALL IMAGEDISTRIBUTION
               IF (MOREPRINTING) WRITE(*,'(I6,2F20.10,F8.2,A,F8.3,2F9.3)') &
               & NITERDONE,ETOTAL/NIMAGE,RMS,AVDEV,'%',SEPARATION,STEPTOT,SQRT(DOT_PRODUCT(V,V)/(NIMAGE*NOPT))

               IF (SQVVGUESS) THEN
                    IF (RMS < SQVVGUESSRMSTOL .OR. NITERDONE > NITERSQVVGUESSMAX) EXIT
               ELSEIF ( RMS<=RMSTOL.AND.NITERDONE>1.AND.NITERDONE>NITERMIN ) THEN
                    EXITSTATUS=1
                    EXIT
               ENDIF

               IF (DEBUG) THEN
                    IF (MOD(NITERDONE,10)==0) CALL DUMPFILES("M")
                    !IF (TWOD) CALL WRITEXYZTMP
               ENDIF
               IF (DUMPNEBXYZ.AND.MOD(NITERDONE,DUMPNEBXYZFREQ)==0) CALL RWG("W",.FALSE.,NITERDONE)
               IF (DUMPNEBPTS.AND.MOD(NITERDONE,DUMPNEBPTSFREQ)==0) CALL SAVEBANDCOORD
               IF (DUMPNEBEOS.AND.MOD(NITERDONE,DUMPNEBEOSFREQ)==0) CALL WRITEPROFILE(NITERDONE)

               ! WE NEED TO DUMP VELOCITY VECTOR
               OPEN(UNIT=333,FILE="VEL",STATUS="REPLACE",FORM="UNFORMATTED")
               WRITE(UNIT=333) V
               CLOSE(333)
          ENDDO
          NITERDONE = NITERDONE - 1

          IF (NITERDONE == NITERMAX) EXITSTATUS=2
          IF (SQVVGUESS) THEN
               NITERDONESAVE = NITERDONE
               INTSTR = WI(NITERDONE)
               PRINT *, 'INITIAL RELAXATION: '//TRIM(INTSTR)//' SQVV STEPS.'
          ENDIF
          IF (DEBUG.AND..NOT.SQVVGUESS) CALL DUMPFILES("E")
          DEALLOCATE(V,DX)

          CONTAINS

          SUBROUTINE QUENCH
               QUENCHING = "J"
               SELECT CASE(QUENCHING)
               CASE ("J") ! QUENCHING BY JONSSON ET AL.
                    DO I=1,NIMAGE
                         PROJ = DOT_PRODUCT(V(NOPT*(I-1)+1:NOPT*I),G(NOPT*(I-1)+1:NOPT*I)) ! NE ZOVISIM PROJ BO G NE NORMOVANY
                         IF (PROJ > 0.0D0) THEN
                              !PRINT '(1X,A,I3,A)', 'FOR IMAGE',I,' VELOCITY WAS ZEROED'
                              V(NOPT*(I-1)+1:NOPT*I) = 0.0D0
                         ELSE
                              GNORM = DOT_PRODUCT(G(NOPT*(I-1)+1:NOPT*I),G(NOPT*(I-1)+1:NOPT*I)) ! NORMA V KVADRATI (DYV *)! 
                              V(NOPT*(I-1)+1:NOPT*I) = -(PROJ/GNORM)*G(NOPT*(I-1)+1:NOPT*I)
                         ENDIF
                    ENDDO
               CASE ("C") ! CREHUET AND FIELD QUENCHING
                    Q = 3 ! ROZMIR KOMIRKY :)
                    DO I=1,NIMAGE*NOPT-(Q-1),Q
                         PROJ = DOT_PRODUCT(V(I:I+Q-1),G(I:I+Q-1))
                         IF (PROJ > 0.0D0) THEN
                              V(I:I+Q-1) = 0.0D0
                         ELSE
                              GNORM = DOT_PRODUCT(G(I:I+Q-1),G(I:I+Q-1))
                              V(I:I+Q-1) = -(PROJ/GNORM)*G(I:I+Q-1)
                         ENDIF
                    ENDDO
               CASE ("G") ! NEIGHBOURS GRADIENTS ARE USED IN PROJECTIONS; THIS IS FOR 2D ONLY
                    Q = 2
                    DO I=1,NIMAGE*NOPT-(2*Q-1),Q
                         PROJ = DOT_PRODUCT(V(I:I+Q-1),G(I:I+Q-1)+G(I-Q:I-1)+G(I+Q:I+2*Q-1))
                         IF (PROJ > 0.0D0) THEN
                              V(I:I+Q-1) = 0.0D0
                         ELSE
                              GNORM = DOT_PRODUCT(G(I:I+Q-1)+G(I-Q:I-1)+G(I+Q:I+2*Q-1),G(I:I+Q-1)+G(I-Q:I-1)+G(I+Q:I+2*Q-1))
                              V(I:I+Q-1) = -(PROJ/GNORM)*(G(I:I+Q-1)+G(I-Q:I-1)+G(I+Q:I+2*Q-1))
                         ENDIF
                    ENDDO
               END SELECT
          END SUBROUTINE QUENCH
     END SUBROUTINE NEBSQVV
     
END MODULE MINIMISER2
