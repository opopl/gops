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
MODULE NEBUTILS
     IMPLICIT NONE
     CONTAINS

     SUBROUTINE IMAGEDISTRIBUTION
          USE NEBDATA
          USE KEYNEB,ONLY: NIMAGE
          IMPLICIT NONE
     
          AVDEVOLD=AVDEV
          
          CALL DISTANCES
          DEVIATION = ABS( 100*( (NIMAGE+1)*DVEC/SEPARATION -1 ) )
          AVDEV     = SUM(DEVIATION)/(NIMAGE+1)
     END SUBROUTINE IMAGEDISTRIBUTION

     SUBROUTINE DISTANCES     !THIS SUBROUTINE IS CALCULATING SEPARATION OF END POINTS AND
          USE KEY,ONLY: BULKT, TWOD
          USE NEBDATA
          USE KEYNEB,ONLY: NIMAGE
          USE COMMONS,ONLY: PARAM1,PARAM2,PARAM3
          IMPLICIT NONE       !DISTANCES BETWEEN INDIVIDUAL IMAGES BOTH
          INTEGER :: J, K

          ! /------------------SEPARATION(ALONG THE MEP)-------------------\
          ! Q--------IMAGE1------IMAGE2---- ...----IMAGENIMAGE-------------FIN
          ! \_DVEC(1)_/ \_DVEC(2)_/                       \_DVEC(NIMAGE+1)_/
     
          IF (BULKT) THEN
             DO J=1,NIMAGE+1
                DVEC(J)=0.0D0
                DO K=1,NATOMS
                   DVEC(J)=DVEC(J)+( XYZ(NOPT*J+3*(K-1)+1) - XYZ(NOPT*(J-1)+3*(K-1)+1) &
  &                    -PARAM1*NINT((XYZ(NOPT*J+3*(K-1)+1) - XYZ(NOPT*(J-1)+3*(K-1)+1))/PARAM1) )**2 &
  &                               +( XYZ(NOPT*J+3*(K-1)+2) - XYZ(NOPT*(J-1)+3*(K-1)+2 ) &
  &                    -PARAM2*NINT((XYZ(NOPT*J+3*(K-1)+2) - XYZ(NOPT*(J-1)+3*(K-1)+2))/PARAM2) )**2 
                   IF (.NOT.TWOD) DVEC(J)=DVEC(J) &
                                  +( XYZ(NOPT*J+3*(K-1)+3) - XYZ(NOPT*(J-1)+3*(K-1)+3 ) &
  &                    -PARAM3*NINT((XYZ(NOPT*J+3*(K-1)+3) - XYZ(NOPT*(J-1)+3*(K-1)+3))/PARAM3) )**2
                ENDDO
             ENDDO
          ELSE
             DO J=1,NIMAGE+1
                  DVEC(J)   = SUM( ( XYZ( NOPT*J+1:NOPT*(J+1) ) - XYZ( NOPT*(J-1)+1:NOPT*J ) )**2)
             ENDDO
          ENDIF
          DVEC(1:NIMAGE+1)=SQRT(DVEC(1:NIMAGE+1))
          SEPARATION=SUM(DVEC(1:NIMAGE+1))

     END SUBROUTINE DISTANCES

     SUBROUTINE INTERNALIMAGEDISTRIBUTION(XINT)
          USE KEY,ONLY: BULKT
          USE NEBDATA
          USE KEYNEB,ONLY: NIMAGE
          IMPLICIT NONE

          DOUBLE PRECISION :: XINT(:),DUMMY
          INTEGER :: J,J1

          DOUBLE PRECISION :: PI=3.141592653589793D0

          AVDEVOLD=AVDEV

          IF (BULKT) THEN
             PRINT '(A)','INTERNALIMAGEDISTRIBUTION> ERROR - THIS ROUTINE SHOULD NOT BE CALLED WITH BULKT TRUE'
             STOP
          ENDIF
          DO J=1,NIMAGE+1
             DVEC(J)=0.0D0
             DO J1=1,NINTS
                DUMMY = XINT( NINTS*J + J1 ) - XINT( NINTS*(J-1) + J1)
                IF (DUMMY.GT.PI) DUMMY=DUMMY-2.0D0*PI
                IF (DUMMY.LT.PI) DUMMY=DUMMY+2.0D0*PI
                DVEC(J)=DVEC(J)+DUMMY*DUMMY
             ENDDO
          ENDDO
          DVEC=SQRT(DVEC)
          SEPARATION=SUM(DVEC)

          DEVIATION = ABS( 100*( (NIMAGE+1)*DVEC/SEPARATION -1 ) )
          AVDEV     = SUM(DEVIATION)/(NIMAGE+1)
     END SUBROUTINE INTERNALIMAGEDISTRIBUTION

     SUBROUTINE MAKEIMAGE(EINITIAL,EFINAL,QQ,FINFIN) ! INTERPOLATES THE BAND
          USE NEBDATA
          USE SPFUNCTS
!          USE KEYNEB,ONLY: NIMAGE
          USE KEY,ONLY: MORPHT, MSTEPS, DEBUG, GREATCIRCLET, GCIMAGE, GCCONV, GCUPDATE, GCSTEPS, INTEPSILON, &
                        BULKT, RBAAT, NTSITES, AMBERT, LOCALPERMDIST, NRBGROUP, RBGROUP, RBNINGROUP, NRBTRIES, NABT,TWOD, RIGIDBODY
          USE INTCOMMONS, ONLY : NNZ, KD, NINTC, DESMINT, DIHINFO, PREVDIH, BACKTCUTOFF, INTERPBACKTCUT, MINBACKTCUT, &
                                 CHICDNEB
          USE COMMONS, ONLY : ZSYM, PARAM1,PARAM2,PARAM3, NRBSITES
          USE MODCHARMM, ONLY : CHRMMT, ICINTERPT
          USE MODAMBER9, ONLY: AMBERICT, AMBICDNEBT, NICTOT
          USE PORFUNCS 
          USE KEYNEB,ONLY: NIMAGE,XYZFILE,RBXYZFILE
          USE KEY, ONLY: FILTH,FILTHSTR  
          IMPLICIT NONE
          INTEGER :: I, J1, ITDONE, INTERVAL, NDONE, J2, J, K, ISTAT, J3, J4, J5, NDUMMY, J6, J7, JWORST2, JWORST3
          DOUBLE PRECISION,ALLOCATABLE :: DELTAX(:), QTN(:,:), PTN(:,:)
          DOUBLE PRECISION DPRAND, SHIFT(NOPT), DUMMY, DLENGTH, DUMMY2, ENERGY, VNEW(NOPT), LRMS, &
   &                       QS(NOPT), QF(NOPT), EINITIAL, EFINAL, COORDS(NOPT), SFRAC
          DOUBLE PRECISION THETA, THETAH, ST, CT, P(3), FCT, XSHIFT, YSHIFT, ZSHIFT, EWORST
          DOUBLE PRECISION, ALLOCATABLE :: XINITIAL(:), XIMAGE(:,:)
                    DOUBLE PRECISION,DIMENSION(:)         :: QQ,FINFIN
          LOGICAL KNOWE, KNOWG, KNOWH, MFLAG, ALIGNCOMMON, ALIGNEDBEFORE(NATOMS)
          COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
          INTEGER NRBSET, NRBTOTAL, NRUNNING, NTRIES, JWORST, ALIGNATOM(NATOMS), NALIGNATOM
          DOUBLE PRECISION CMXS, CMYS, CMZS, CMXF, CMYF, CMZF, LSTART(3*NATOMS), LFINISH(3*NATOMS), &
  &                        RBANGLE, E1, E2, LTEMP(3*NATOMS), LEIMAGE(NIMAGE), LEIMAGE2(NIMAGE), LEIMAGE3(NIMAGE), &
  &                        LPRED(3*NATOMS), THETA1, THETA2, LX(3), LV(3), DBEST, TBEST, LVBEST(3), LGDUMMY(3*NATOMS), RMSDUMMY, &
  &                        QS2(3*NATOMS), QF2(3*NATOMS), D, EWORST2, EWORST3

          DOUBLE PRECISION :: STXYZ(3*NTSITES)  
          CHARACTER(LEN=80) :: FILENAME,FILENAME2,DUMMYS,DUMMYS2
           
          ! FOR INTERNALS
          DOUBLE PRECISION :: DELTACART(3*NATOMS)
          LOGICAL :: INTPTEST, FAILED
          ! ALIGNMENT STUFF
          DOUBLE PRECISION :: DISTF, DIST, DIST2, RMAT(3,3)
          CHARACTER(LEN=5) :: ZSYMSAVE
          COMMON /SYS/ ZSYMSAVE

          IF (DESMINT) THEN
             BACKTCUTOFF = INTERPBACKTCUT
             INTPTEST = .FALSE.
          ENDIF

          IF (MORPHT) THEN
             IF (DESMINT) THEN
                PRINT*, 'ERROR! MORPHT NOT IMPLEMENTED WITH DESM.'
                STOP
             ENDIF

             QS(1:3*NATOMS)=XYZ(1:NOPT)
             QF(1:3*NATOMS)=XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+1)+NOPT)
             KNOWE=.FALSE.;  KNOWG=.FALSE.; KNOWH=.FALSE.
          
             CALL MORPH(MSTEPS,QS,QF,ENERGY,VNEW,MFLAG,LRMS,ITDONE,.TRUE.)
             IF (ITDONE.LT.NIMAGE) THEN
                IF (MFLAG) THEN
!
!  UNFORTUNATELY MANY STATEMENTS IN THE NEB
!  ROUTINES DO NOT SPECIFY THE ARRAY BOUNDS ON ARRAY OPERATIONS, WHICH THEN DEFAULT TO THE
!  DECLARED ARRAY SIZE AND THEREFORE GO WRONG. HENCE IF WE REDUCE THE NUMBER OF IMAGES HERE
!  THEN WE HAVE TO DEALLOCATE AND REALLOCATE:
!                  
                   PRINT '(A,I6)',' MAKEIMAGE> FEWER MORPH POINTS THAN IMAGES - REDUCING IMAGES TO ',ITDONE
                   NIMAGE=ITDONE

                   DEALLOCATE(XYZ,GGG,SSS,EEE,RRR,TANVEC,DVEC,DEVIATION,STEPIMAGE,TRUEGRAD)
                   ALLOCATE(XYZ(NOPT*(NIMAGE+2)),GGG(NOPT*(NIMAGE+2)),SSS(NOPT*(NIMAGE+2)),EEE(NIMAGE+2), &
  &   RRR(NIMAGE+2),TANVEC(NOPT,NIMAGE),DVEC(NIMAGE+1),DEVIATION(NIMAGE+1),STEPIMAGE(NIMAGE),TRUEGRAD(NOPT*(NIMAGE+2)))
                   X         => XYZ(NOPT+1:NOPT*(NIMAGE+1))
                   EIMAGE    => EEE(2:NIMAGE+1)
                   G         => GGG(NOPT+1:NOPT*(NIMAGE+1))
                   GSPR      => SSS(NOPT+1:NOPT*(NIMAGE+1))
                   RMSFIMAGE => RRR(2:NIMAGE+1)
                   EEE = 0.0D0
                   EEE(1)=EINITIAL
                   EEE(NIMAGE+2)=EFINAL
                   XYZ(:NOPT)=QQ
                   XYZ(NOPT*(NIMAGE+1)+1:)=FINFIN

                   INTERVAL=ITDONE/NIMAGE
                   NDONE=0
                   OPEN(UNIT=991,FILE='MORPH.POINTS',STATUS='OLD')
                   DO J1=2,NIMAGE+1
                      DO J2=1,INTERVAL
                         NDONE=NDONE+1
                         READ(991,*) XYZ(NOPT*(J1-1)+1:NOPT*(J1-1)+NOPT) 
                      ENDDO
                      IF (DEBUG) PRINT '(2(A,I6))',' MAKEIMAGE> IMAGE ',J1,' READ FROM MORPH.POINTS FRAME ',NDONE
                   ENDDO
                   CLOSE(991)
                   XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+1)+NOPT)=QF(1:NOPT) ! TO ALIGN THE FINAL IMAGE
                ELSE
                   PRINT '(A)',' MAKEIMAGE> FEWER MORPH POINTS THAN IMAGES - REVERT TO LINEAR INTERPOLATION'
                   ALLOCATE(DELTAX(NOPT))
                   DELTAX(1:NOPT) = ( XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+2) ) - XYZ(1:NOPT) )/(NIMAGE+1)
                   DO I=2,NIMAGE+1
                        XYZ(NOPT*(I-1)+1:NOPT*I) = XYZ(1:NOPT) + DELTAX*(I-1)
                   ENDDO
                   DEALLOCATE(DELTAX)
                ENDIF
             ELSE
                INTERVAL=ITDONE/NIMAGE
                NDONE=0
                OPEN(UNIT=991,FILE='MORPH.POINTS',STATUS='OLD')
                DO J1=2,NIMAGE+1
                   DO J2=1,INTERVAL
                      NDONE=NDONE+1
                      READ(991,*) XYZ(NOPT*(J1-1)+1:NOPT*(J1-1)+NOPT) 
                   ENDDO
                   IF (DEBUG) PRINT '(2(A,I6))',' MAKEIMAGE> IMAGE ',J1,' READ FROM MORPH.POINTS FRAME ',NDONE
                ENDDO
                CLOSE(991)
                XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+1)+NOPT)=QF(1:NOPT) ! TO ALIGN THE FINAL IMAGE
             ENDIF
             MORPHT=.FALSE. ! DJW
          ELSEIF (GREATCIRCLET) THEN 
             IF (NOPT.NE.3*NATOMS) THEN
                PRINT '(A)','NNUTILS> ERROR - NOPT NEEDS TO BE 3*NATOMS TO USE GREATCIRCLE INTERPOLATION'
                STOP
             ENDIF
             ALLOCATE(XINITIAL(3*NATOMS+1),XIMAGE(NIMAGE,3*NATOMS))
             XINITIAL(1:3*NATOMS+1)= 0.0D0
             XINITIAL(3*NATOMS)= 1.0D0
             XINITIAL(3*NATOMS+1)= 1.0D0
             GCIMAGE=NIMAGE
             QS(1:3*NATOMS)=XYZ(1:NOPT)
             QF(1:3*NATOMS)=XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+1)+NOPT)
             CALL GCLBFGS(QS,QF,XIMAGE,3*NATOMS+1,GCUPDATE,XINITIAL,.FALSE.,GCCONV,MFLAG,ENERGY,RMS,GCSTEPS,.TRUE.,ITDONE,.TRUE.)
             DO I=2,NIMAGE+1
                  XYZ(NOPT*(I-1)+1:NOPT*I) = XIMAGE(I-1,1:NOPT)
             ENDDO
             DEALLOCATE(XINITIAL,XIMAGE)

!     DC430 >
          ELSEIF (RBAAT) THEN

             ALLOCATE(DELTAX(NOPT))
             ALLOCATE(QTN(NIMAGE+2,4))
             ALLOCATE(PTN(NIMAGE+2,4))

             DELTAX(:) = (XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+2)) - XYZ(1:NOPT))/(NIMAGE+1)
            
!             PRINT *, XYZ(9), XYZ(NOPT*(NIMAGE+1)+9) 
!             PRINT *, 9, DELTAX(9) 
             DO I = 2, NIMAGE+1
                J = NOPT*(I-1)
                XYZ(J+1:J+NOPT/2) = XYZ(1:NOPT/2) + DELTAX(1:NOPT/2)*(I-1)
!                PRINT *, I, XYZ(J+9), XYZ(J+9) - XYZ(J+9-NOPT)
             ENDDO

             DO J1 = 1, NATOMS/2

!     CONVERT FROM AA -> QUATERNION FOR THE INITIAL AND FINAL FRAMES

                DO K = 1, (NIMAGE+2), (NIMAGE+1)

                   J = (K-1)*NOPT + NOPT/2 + 3*(J1-1)
                   P(:)   = XYZ(J+1:J+3)
                   THETA  = DSQRT(DOT_PRODUCT(P(:),P(:)))
                   IF (THETA == 0.D0) THEN
                      QTN(K,1)   = 1.D0
                      QTN(K,2:4) = 0.D0
                   ELSE
                      THETAH     = 0.5D0*THETA
                      ST         = SIN(THETAH)
                      QTN(K,1)   = COS(THETAH)
                      QTN(K,2:4) = P(:)*ST/THETA
                   ENDIF
                ENDDO

!     QUATERNION INTERPOLATION: ISLERP
!     NOW THETA = \ALPHA (THE ANGLE BETWEEN THE TWO QUATERNIONS)

                CT       = DOT_PRODUCT(QTN(1,:),QTN(NIMAGE+2,:))
                IF (CT < 0.D0) THEN
                   CT =-CT
                   QTN(NIMAGE+2,:) =-QTN(NIMAGE+2,:)
                ENDIF
                THETA    = ACOS(CT)
                ST       = SIN(THETA)
!     INCREMENTAL APPROACH: TANGENT QUATERNION
                PTN(1,:) = (QTN(NIMAGE+2,:) - CT*QTN(1,:))/ST
!     NOW THETA = \BETA = \ALPHA/(NIMAGE+1)
                THETA    = THETA/(NIMAGE+1)
                ST       = SIN(THETA)
                CT       = COS(THETA)

                DO I = 2, NIMAGE+1

                   QTN(I,:) = CT*QTN(I-1,:) + ST*PTN(I-1,:)
                   PTN(I,:) = CT*PTN(I-1,:) - ST*QTN(I-1,:)

!     CONVERT FROM QUATERNION -> AA
                   THETA  = 2.D0*ACOS(QTN(I,1))
                   J = NOPT*(I-1) + NOPT/2 + 3*(J1-1)
                   IF (THETA == 0.D0) THEN
                      XYZ(J+1:J+3) = 0.D0
                   ELSE
                      FCT    = DSQRT(DOT_PRODUCT(QTN(I,2:4),QTN(I,2:4)))
                      XYZ(J+1:J+3) = THETA*QTN(I,2:4)/FCT
                   ENDIF
                  
                ENDDO
  
             ENDDO
 
             DEALLOCATE(DELTAX)
             DEALLOCATE(QTN)
             DEALLOCATE(PTN)

          ELSE
             ALLOCATE(DELTAX(NOPT))

             IF (BULKT) THEN
                DO K=1,NATOMS
                   DELTAX(3*(K-1)+1)=XYZ(NOPT*(NIMAGE+1)+3*(K-1)+1) - XYZ(3*(K-1)+1) &
  &                    -PARAM1*NINT((XYZ(NOPT*(NIMAGE+1)+3*(K-1)+1) - XYZ(3*(K-1)+1))/PARAM1) 
                   DELTAX(3*(K-1)+2)=XYZ(NOPT*(NIMAGE+1)+3*(K-1)+2) - XYZ(3*(K-1)+2) &
  &                    -PARAM2*NINT((XYZ(NOPT*(NIMAGE+1)+3*(K-1)+2) - XYZ(3*(K-1)+2))/PARAM2)
                   IF (.NOT.TWOD) DELTAX(3*(K-1)+3)=XYZ(NOPT*(NIMAGE+1)+3*(K-1)+3) - XYZ(3*(K-1)+3) &
  &                    -PARAM3*NINT((XYZ(NOPT*(NIMAGE+1)+3*(K-1)+3) - XYZ(3*(K-1)+3))/PARAM3)
                ENDDO
                DELTAX(1:NOPT)=DELTAX(1:NOPT)/(NIMAGE+1)
             ELSE
                DELTAX(1:NOPT) = ( XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+2) ) - XYZ(1:NOPT) )/(NIMAGE+1)
             ENDIF
   
! BS360: INTERPOLATION USING DIHEDRALS.
             IF (CHRMMT.AND.CHICDNEB.AND.ICINTERPT) THEN
                QS(1:NOPT)=XYZ(1:NOPT)         
                QF(1:NOPT)=XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+2))
             ENDIF
! END BS360
! MSB50: INTERPOLATION USING DIHEDRALS FOR AMBER
             IF ((AMBERT.OR.NABT).AND.AMBICDNEBT.AND.AMBERICT) THEN
                QS(1:NOPT)=XYZ(1:NOPT)
                QF(1:NOPT)=XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+2))
                PRINT*, "MSB50 NOPT", NOPT
!                DO K=1, 3*NATOMS
!                    PRINT '(6F9.3)', QF(6*(K-1)+1:6*K)
!                ENDDO
             ENDIF
! END MSB50
             DO I=2,NIMAGE+1
                  XYZ(NOPT*(I-1)+1:NOPT*I) = XYZ(1:NOPT) + DELTAX*(I-1)
       
! BS360
                  IF (CHRMMT.AND.CHICDNEB.AND.ICINTERPT) THEN
                     COORDS(1:NOPT)=XYZ(NOPT*(I-1)+1:NOPT*I)
                     SFRAC=1.D0-(I-1.D0)/(NIMAGE+1.D0)                  
!                     WRITE(6,*) 'I,NIMAGE+1,SFRAC= ',I,NIMAGE+1,SFRAC
                     CALL ICINTERPOL(COORDS,QS,QF,SFRAC)
                     XYZ(NOPT*(I-1)+1:NOPT*I)=COORDS(1:NOPT)        
!                     WRITE(100+I,'(3F10.3)') XYZ(NOPT*(I-1)+1:NOPT*I)
!                     CALL FLUSH(100+I)
                     CALL FLUSH(100+I,ISTAT)
                  ENDIF
! END BS360

! MSB50
                  IF ((AMBERT.OR.NABT).AND.AMBICDNEBT.AND.AMBERICT) THEN
                     COORDS(1:NOPT)=XYZ(NOPT*(I-1)+1:NOPT*I)
                     SFRAC=1.D0-(I-1.D0)/(NIMAGE+1.D0)
                     IF (.NOT.ALLOCATED(NICTOT)) THEN
                         CALL SETDIHEAM() 
                     ENDIF
                     CALL TAKESTEPAMDIHED(COORDS, QS, QF,SFRAC)
                     WRITE(6,*) 'I,NIMAGE+1,SFRAC= ',I,NIMAGE+1,SFRAC
                     XYZ(NOPT*(I-1)+1:NOPT*I)=COORDS(1:NOPT)
                     !DO K=1,(NOPT/3)
                     !  PRINT*, XYZ(3*(K-1)+1:3*K) 
                     !ENDDO
!                    WRITE(100+I,'(3F10.3)') XYZ(NOPT*(I-1)+1:NOPT*I)
!                     CALL FLUSH(100+I)
                  ENDIF
! END MSB50

                  IF (DESMINT) THEN
                     PREVDIH => DIHINFO(I,:)
                     CALL TRANSBACKDELTA(DELTAX,DELTACART,XYZCART(3*NATOMS*(I-2)+1:3*NATOMS*(I-1)),&
                          & NINTC,3*NATOMS,NNZ,KD,FAILED,INTPTEST,INTEPSILON)
                     XYZCART(3*NATOMS*(I-1)+1:3*NATOMS*I) = XYZCART(3*NATOMS*(I-2)+1:3*NATOMS*(I-1)) + DELTACART(1:3*NATOMS)
                     CALL NEWMINDIST(XYZCART(1:3*NATOMS),XYZCART(3*NATOMS*(I-1)+1:3*NATOMS*I),&
                          & NATOMS,DIST,.FALSE.,.FALSE.,ZSYM(1),.FALSE.,.FALSE.,.FALSE.,RMAT)
                  ENDIF
             ENDDO
             
             IF (DESMINT) BACKTCUTOFF = MINBACKTCUT

             DEALLOCATE(DELTAX)
          ENDIF
!
!  NOW IF WE HAVE SOME RIGID BODIES, INTERPOLATE THEM AS RIGID BODIES!
!  HOWEVER, ALSO CHECK WHETHER THE ALTERNATIVE INTERPOLATION BASED ON
!  REGULAR MINPERMDIST FOR THE WHOLE MOLECULE GIVES A BETTER INITIAL GUESS.
!  IF SO, WE CAN REPLACE BLOCKS OF ATOMS WITH THIS INTERPOLATION, BUT
!  WE MUST USE THE SAME LOCAL PERMUTATIONAL ISOMER FOR EACH IMAGE.
!
     IF (LOCALPERMDIST) THEN
        ALIGNEDBEFORE(1:NATOMS)=.FALSE.
        IF (NRBGROUP.GT.0) THEN
           NRBTOTAL=0
           DO J1=1,NRBGROUP
              NRBTOTAL=NRBTOTAL+RBNINGROUP(J1)
           ENDDO
           NTRIES=0
11         CONTINUE
!
!  IDENTIFY THE CURRENT HIGHEST IMAGE. WE WILL USE THIS TO CHOOSE BETWEEN
!  + AND - ROTATIONS ABOUT THE CHOSEN RB AXIS FOR EACH IMAGE.
!
           EWORST=-HUGE(1.0D0)
           DO J1=1,NIMAGE
              COORDS(1:3*NATOMS)=XYZ(NOPT*J1+1:NOPT*J1+3*NATOMS)
              CALL POTENTIAL(COORDS,E1,LGDUMMY,.FALSE.,.FALSE.,RMSDUMMY,.FALSE.,.FALSE.)
              LEIMAGE(J1)=E1
              IF (E1.GT.EWORST) THEN
                 JWORST=J1
                 EWORST=E1
              ENDIF
!             IF (DEBUG) PRINT '(A,I6,A,G20.10)',' NNUTILS> IMAGE ',J1,' ENERGY=',E1
           ENDDO
           IF (DEBUG) PRINT '(A,I6,A,G20.10)',' NNUTILS> HIGHEST IMAGE IS ',JWORST,' ENERGY=',EWORST
           NTRIES=NTRIES+1
!          PRINT '(A,I6)',' MAKEIMAGE> NUMBER OF ENTRIES IN RBGROUP=',NRBTOTAL
!          PRINT '(22I6)',RBGROUP(1:NRBTOTAL)
!
!  RIGID BODY INTERPOLATIONS USE LSTART AND LFINISH, WHICH ARE ASSIGNED FROM QS AND QF.
!  WE REPLACE BLOCKS OF ATOMS IN THE XYZ ARRAY CORRESPONDING TO THESE ATOMS IF THE LOCAL
!  RIGID BODY INTERPOLATION LOWERS THE ENERGY.
!
           QS(1:3*NATOMS)=XYZ(1:NOPT)
           QF(1:3*NATOMS)=XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+1)+NOPT)
!
!  NOW CHANGE THE WHOLE XYZ ARRAY TO BE BASED ON THE PERMUTATIONAL ISOMER THAT
!  MINIMISES THE OVERALL DISTANCE. THE DIFFERENCE IS IN THE INITIAL STRAIGHT LINE
!  GUESS, WHICH MAY BE FOR DIFFERENT PERMUTATIONAL ISOMERS.
!  USE THE PERMUTATION THAT GIVES THE LOWEST VALUE FOR THE HIGHEST IMAGE.
!
           QS2(1:3*NATOMS)=XYZ(1:NOPT)
           QF2(1:3*NATOMS)=XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+1)+NOPT)
           CALL MINPERMDIST(QS2,QF2,NATOMS, &
  &                         DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
           EWORST2=-HUGE(1.0D0)
           DO J1=1,NIMAGE
              COORDS(1:3*NATOMS)=(J1*QF2(1:3*NATOMS)+(NIMAGE+1-J1)*QS2(1:3*NATOMS))/(NIMAGE+1)
              CALL POTENTIAL(COORDS,E1,LGDUMMY,.FALSE.,.FALSE.,RMSDUMMY,.FALSE.,.FALSE.)
              LEIMAGE2(J1)=E1
              IF (E1.GT.EWORST2) THEN
                 JWORST2=J1
                 EWORST2=E1
              ENDIF
!             IF (DEBUG) PRINT '(A,I6,A,G20.10)',' NNUTILS> ALTERNATIVE INTERPOLATION, IMAGE ',J1,' ENERGY=',E1
              IF (EWORST2.GT.EWORST) EXIT ! NO POINT CONTINUING, WE WILL USE THE OTHER INITIAL GUESS.
           ENDDO
           IF (EWORST2.LT.EWORST) THEN
              IF (DEBUG) PRINT '(A)',' NNUTILS> ADOPTING ALTERNATIVE INTERPOLATION FOR INITIAL GUESS'
              JWORST=JWORST2
              EWORST=EWORST2 
              LEIMAGE(1:NIMAGE)=LEIMAGE2(1:NIMAGE)
              IF (DEBUG) PRINT '(A,I6,A,G20.10)',' NNUTILS> HIGHEST IMAGE IS ',JWORST,' ENERGY=',EWORST
              XYZ(1:3*NATOMS)=QS2(1:3*NATOMS)
              XYZ((NIMAGE+1)*NOPT+1:(NIMAGE+1)*NOPT+3*NATOMS)=QF2(1:3*NATOMS)
              DO J1=1,NIMAGE
                 XYZ(J1*NOPT+1:J1*NOPT+3*NATOMS)=(J1*QF2(1:3*NATOMS)+(NIMAGE+1-J1)*QS2(1:3*NATOMS))/(NIMAGE+1)
              ENDDO
           ENDIF
!
! INITIAL GUESS HAS BEEN SELECTED
!
           NRUNNING=0
           GROUPLOOP: DO J1=1,NRBGROUP
              CMXS=0.0D0; CMYS=0.0D0; CMZS=0.0D0; CMXF=0.0D0; CMYF=0.0D0; CMZF=0.0D0
              DO J2=1,RBNINGROUP(J1)
                 J3=RBGROUP(J2+NRUNNING)
                  LSTART(3*(J2-1)+1:3*(J2-1)+3)=QS(3*(J3-1)+1:3*(J3-1)+3)
                 LFINISH(3*(J2-1)+1:3*(J2-1)+3)=QF(3*(J3-1)+1:3*(J3-1)+3)
                 CMXS=CMXS+LSTART(3*(J2-1)+1)
                 CMYS=CMYS+LSTART(3*(J2-1)+2)
                 CMZS=CMZS+LSTART(3*(J2-1)+3)
                 CMXF=CMXF+LFINISH(3*(J2-1)+1)
                 CMYF=CMYF+LFINISH(3*(J2-1)+2)
                 CMZF=CMZF+LFINISH(3*(J2-1)+3)
!                PRINT '(3(A,I6))',' MAKEIMAGE> GROUP ',J1,' CONTAINS ATOM ',J3,' TOTAL IN GROUP=',RBNINGROUP(J1)
              ENDDO
              NRBSET=RBNINGROUP(J1)
              NRUNNING=NRUNNING+NRBSET
              CMXS=CMXS/NRBSET; CMYS=CMYS/NRBSET; CMZS=CMZS/NRBSET;
              CMXF=CMXF/NRBSET; CMYF=CMYF/NRBSET; CMZF=CMZF/NRBSET;
              DO J2=1,NRBSET
                 LSTART(3*(J2-1)+1)=LSTART(3*(J2-1)+1)-CMXS
                 LSTART(3*(J2-1)+2)=LSTART(3*(J2-1)+2)-CMYS
                 LSTART(3*(J2-1)+3)=LSTART(3*(J2-1)+3)-CMZS
                 LFINISH(3*(J2-1)+1)=LFINISH(3*(J2-1)+1)-CMXF
                 LFINISH(3*(J2-1)+2)=LFINISH(3*(J2-1)+2)-CMYF
                 LFINISH(3*(J2-1)+3)=LFINISH(3*(J2-1)+3)-CMZF
              ENDDO
              IF (DEBUG) THEN
                 PRINT '(A,I8,A,I8)',' MAKEIMAGE> NUMBER OF ATOMS IN RIGID BODY GROUP ',J1,' IS ',NRBSET
!                PRINT '(A,3F15.5)',' MAKEIMAGE> START  CENTRE OF MASS: ',CMXS,CMYS,CMZS
!                PRINT '(A,3F15.5)',' MAKEIMAGE> START  COORDINATES:'
!                PRINT '(3F20.10)',LSTART(1:3*NRBSET)
!                PRINT '(A,3F15.5)',' MAKEIMAGE> FINISH CENTRE OF MASS: ',CMXF,CMYF,CMZF
!                PRINT '(A,3F15.5)',' MAKEIMAGE> FINISH COORDINATES:'
!                PRINT '(3F20.10)',LFINISH(1:3*NRBSET)
              ENDIF
              DBEST=1.0D100
              TBEST=0.0D0
              DO J2=1,NRBSET
                 DO J3=J2+1,NRBSET
                    LV(1)= (LSTART(3*(J2-1)+2)-LFINISH(3*(J2-1)+2))*(LSTART(3*(J3-1)+3)-LFINISH(3*(J3-1)+3)) &
  &                     -(LSTART(3*(J2-1)+3)-LFINISH(3*(J2-1)+3))*(LSTART(3*(J3-1)+2)-LFINISH(3*(J3-1)+2))
                    LV(2)=-(LSTART(3*(J2-1)+1)-LFINISH(3*(J2-1)+1))*(LSTART(3*(J3-1)+3)-LFINISH(3*(J3-1)+3)) &
  &                     +(LSTART(3*(J2-1)+3)-LFINISH(3*(J2-1)+3))*(LSTART(3*(J3-1)+1)-LFINISH(3*(J3-1)+1))
                    LV(3)= (LSTART(3*(J2-1)+1)-LFINISH(3*(J2-1)+1))*(LSTART(3*(J3-1)+2)-LFINISH(3*(J3-1)+2)) &
  &                     -(LSTART(3*(J2-1)+2)-LFINISH(3*(J2-1)+2))*(LSTART(3*(J3-1)+1)-LFINISH(3*(J3-1)+1)) 
                    DUMMY=LV(1)**2+LV(2)**2+LV(3)**2
                    IF (DUMMY.EQ.0.0D0) CYCLE
                    DUMMY=SQRT(DUMMY)
                    LV(1)=LV(1)/DUMMY; LV(2)=LV(2)/DUMMY; LV(3)=LV(3)/DUMMY
!
! CHECK RESULT BY BACK-ROTATION AND FIND THE BEST ANGLE-AXIS THAT GIVES THE MINIMUM ERROR.
!
! DOT PRODUCT IN DUMMY SHOULD BE THE SAME FOR THE ATOM IN START AND FINISH.
! HOWEVER, WE ARE DEALING WITH APPROXIMATE REIGID BODIES, SO THE SYMMETRY COULD BE BROKEN
! WITH A LIKELY LOSS OF ACCURACY.
!
                    DUMMY= LSTART(3*(J2-1)+1)*LV(1)+LSTART(3*(J2-1)+2)*LV(2)+LSTART(3*(J2-1)+3)*LV(3)
                    DUMMY=LSTART(3*(J2-1)+1)**2+LSTART(3*(J2-1)+2)**2+LSTART(3*(J2-1)+3)**2-DUMMY**2
                    DUMMY2=(LSTART(3*(J2-1)+1)-LFINISH(3*(J2-1)+1))**2 &
  &                       +(LSTART(3*(J2-1)+2)-LFINISH(3*(J2-1)+2))**2 &
  &                       +(LSTART(3*(J2-1)+3)-LFINISH(3*(J2-1)+3))**2
                    DUMMY2=SQRT(DUMMY2)
                    IF (ABS(DUMMY2/(2.0D0*SQRT(MAX(DUMMY,1.0D-30)))).GT.1.0D0) CYCLE
                    THETA1=2.0D0*ASIN(DUMMY2/(2.0D0*SQRT(MAX(DUMMY,1.0D-30))))

                    DUMMY=LSTART(3*(J3-1)+1)*LV(1)+LSTART(3*(J3-1)+2)*LV(2)+LSTART(3*(J3-1)+3)*LV(3)
                    DUMMY=LSTART(3*(J3-1)+1)**2+LSTART(3*(J3-1)+2)**2+LSTART(3*(J3-1)+3)**2-DUMMY**2
                    DUMMY2=(LSTART(3*(J3-1)+1)-LFINISH(3*(J3-1)+1))**2 &
  &                       +(LSTART(3*(J3-1)+2)-LFINISH(3*(J3-1)+2))**2 &
  &                       +(LSTART(3*(J3-1)+3)-LFINISH(3*(J3-1)+3))**2
                    DUMMY2=SQRT(DUMMY2)
                    IF (ABS(DUMMY2/(2.0D0*SQRT(MAX(DUMMY,1.0D-30)))).GT.1.0D0) CYCLE
                    THETA2=2.0D0*ASIN(DUMMY2/(2.0D0*SQRT(MAX(DUMMY,1.0D-30))))
!
! TRY THETA1
!
                    DUMMY=0.0D0
                    DO J4=1,NRBSET
                       LX(1:3)=LSTART(3*(J4-1)+1:3*(J4-1)+3)
                       CALL ROTATEANGLAXIS(LX,THETA1,LV)
                       DUMMY=DUMMY+(LX(1)-LFINISH(3*(J4-1)+1))**2+(LX(2)-LFINISH(3*(J4-1)+2))**2+(LX(3)-LFINISH(3*(J4-1)+3))**2
                    ENDDO
                    IF (DUMMY.LT.DBEST) THEN
                       DBEST=DUMMY
                       TBEST=THETA1
                       LVBEST(1:3)=LV(1:3)
                    ENDIF
!                   IF (DEBUG) PRINT '(A,2F20.10,A,3F15.5)',' MAKEIMAGE> THETA1, RESIDUAL: ', &
! &                            THETA1,SQRT(DUMMY),' AXIS: ',LV(1:3)  ! DJW
!
! TRY THETA2
!
                    DUMMY=0.0D0
                    DO J4=1,NRBSET
                       LX(1:3)=LSTART(3*(J4-1)+1:3*(J4-1)+3)
                       CALL ROTATEANGLAXIS(LX,THETA2,LV)
                       DUMMY=DUMMY+(LX(1)-LFINISH(3*(J4-1)+1))**2+(LX(2)-LFINISH(3*(J4-1)+2))**2+(LX(3)-LFINISH(3*(J4-1)+3))**2
                    ENDDO
                    IF (DUMMY.LT.DBEST) THEN
                       DBEST=DUMMY
                       TBEST=THETA2
                       LVBEST(1:3)=LV(1:3)
                    ENDIF
!                   IF (DEBUG) PRINT '(A,2F20.10,A,3F15.5)',' MAKEIMAGE> THETA2, RESIDUAL: ', &
! &                            THETA2,SQRT(DUMMY),' AXIS: ',LV(1:3)  ! DJW
                    THETA1=-THETA1; THETA2=-THETA2
!
! TRY -THETA1
!
                    DUMMY=0.0D0
                    DO J4=1,NRBSET
                       LX(1:3)=LSTART(3*(J4-1)+1:3*(J4-1)+3)
                       CALL ROTATEANGLAXIS(LX,THETA1,LV)
                       DUMMY=DUMMY+(LX(1)-LFINISH(3*(J4-1)+1))**2+(LX(2)-LFINISH(3*(J4-1)+2))**2+(LX(3)-LFINISH(3*(J4-1)+3))**2
                    ENDDO
                    IF (DUMMY.LT.DBEST) THEN
                       DBEST=DUMMY
                       TBEST=THETA1
                       LVBEST(1:3)=LV(1:3)
                    ENDIF
!                   IF (DEBUG) PRINT '(A,2F20.10,A,3F15.5)',' MAKEIMAGE> THETA1, RESIDUAL: ', &
! &                            THETA1,SQRT(DUMMY),' AXIS: ',LV(1:3)  ! DJW
!
! TRY -THETA2
!
                    DUMMY=0.0D0
                    DO J4=1,NRBSET
                       LX(1:3)=LSTART(3*(J4-1)+1:3*(J4-1)+3)
                       CALL ROTATEANGLAXIS(LX,THETA2,LV)
                       DUMMY=DUMMY+(LX(1)-LFINISH(3*(J4-1)+1))**2+(LX(2)-LFINISH(3*(J4-1)+2))**2+(LX(3)-LFINISH(3*(J4-1)+3))**2
                    ENDDO
                    IF (DUMMY.LT.DBEST) THEN
                       DBEST=DUMMY
                       TBEST=THETA2
                       LVBEST(1:3)=LV(1:3)
                    ENDIF
!                   IF (DEBUG) PRINT '(A,2F20.10,A,3F15.5)',' MAKEIMAGE> THETA2, RESIDUAL: ', &
! &                            THETA2,SQRT(DUMMY),' AXIS: ',LV(1:3)  ! DJW
                 ENDDO
              ENDDO
!
! IF THE ANGLE IS SMALL (RADIANS) THEN PRESUMABLY A STRAIGHT LINE INTERPOLATION SHOULD BE OK.
! TBEST IS ALSO SET TO ZERO IF THE ATTEMPTED RIGID BODY ROTATION WAS NOT ACCURATE ENOUGH.
!
              IF (ABS(TBEST).LT.0.1D0) CYCLE GROUPLOOP

              IF (DEBUG) THEN
                 PRINT '(2(A,G20.10))',' MAKEIMAGE> BEST RESIDUAL IS ',SQRT(DBEST),' FOR ANGLE ',TBEST
!                PRINT '(A,3F15.5)',' MAKEIMAGE> ROTATED COORDINATES:'
!                DO J4=1,NRBSET
!                   LX(1:3)=LSTART(3*(J4-1)+1:3*(J4-1)+3)
!                   CALL ROTATEANGLAXIS(LX,TBEST,LVBEST)
!                   DUMMY=DUMMY+(LX(1)-LFINISH(3*(J4-1)+1))**2+(LX(2)-LFINISH(3*(J4-1)+2))**2+(LX(3)-LFINISH(3*(J4-1)+3))**2
!                   PRINT '(3F20.10)',LX(1:3)
!                ENDDO
              ENDIF

!
!  INTERPOLATE THIS GROUP USING THE BEST ROTATION. 
!  WE CAN TRY TBEST/2 AND (2*PI-TBEST)/2 TO TRY AND GUESS THE BEST PATH.
!
!  IF THE GROUP CONTAINS AN ATOM THAT IS IN A PREVIOUS RIGID BODY THEN TRANSLATE THE
!  GROUP SO THAT THIS ATOM LINES UP WITH ITS PREVIOUS IMAGES. OTHERWISE SHIFT THE
!  CENTRE OF MASS SEQUENTIALLY BETWEEN THE VALUES FOR THE STARTING AND FINISHING 
!  GEOMETRIES. IT IS POSSIBLE THAT THE CURRENT GROUP WILL CONTAIN MORE THAN ONE ATOM
!  IN COMMON WITH PREVIOUS GROUPS. IN THIS CASE WE USE THE LAST ATOM.
!
              ALIGNCOMMON=.FALSE.
              NALIGNATOM=0
              ALIGNATOM(1:NATOMS)=-1
              DO J4=1,NRUNNING-NRBSET
                 IF (.NOT.ALIGNEDBEFORE(RBGROUP(J4))) CYCLE
                 DO J5=1,NRBSET
                    IF (RBGROUP(NRUNNING-NRBSET+J5).EQ.RBGROUP(J4)) THEN
                       ALIGNCOMMON=.TRUE.
                       NALIGNATOM=NALIGNATOM+1
                       ALIGNATOM(NALIGNATOM)=J5
                       IF (DEBUG) PRINT '(A,I6,A)',' MAKEIMAGE> ATOM ',RBGROUP(NRUNNING-NRBSET+J5), &
  &                                                ' WAS IN A PREVIOUSLY ALIGNED RIGID GROUP'
                    ENDIF
                 ENDDO
              ENDDO
              IF (DEBUG.AND.ALIGNCOMMON) PRINT '(A,I6,A)',' MAKEIMAGE> TRANSLATING GROUP TO ALIGN EACH SHARED ATOM IN TURN '

              ALIGNLOOP: DO J7=1,MAX(NALIGNATOM,1)
                 IF (DEBUG.AND.ALIGNCOMMON) PRINT '(A,I6,A)',' MAKEIMAGE> ALIGNING ON ATOM ',RBGROUP(NRUNNING-NRBSET+ALIGNATOM(J7))

                 EWORST2=-HUGE(1.0D0)
                 DO J4=1,NIMAGE
                    COORDS(1:NOPT)=XYZ(NOPT*J4+1:NOPT*J4+NOPT) 

                    IF (J4.LE.NIMAGE/2) THEN 
                       THETA=J4*TBEST/(NIMAGE+1)
                    ELSE
                       THETA=-(NIMAGE-J4+1)*TBEST/(NIMAGE+1)
                    ENDIF
                    IF (ALIGNCOMMON) THEN
                       IF (J4.LE.NIMAGE/2) THEN
                          LX(1:3)=LSTART(3*(ALIGNATOM(J7)-1)+1:3*(ALIGNATOM(J7)-1)+3)
                       ELSE
                          LX(1:3)=LFINISH(3*(ALIGNATOM(J7)-1)+1:3*(ALIGNATOM(J7)-1)+3)
                       ENDIF
                       CALL ROTATEANGLAXIS(LX,THETA,LVBEST)
                       XSHIFT=XYZ(NOPT*J4+3*(RBGROUP(NRUNNING-NRBSET+ALIGNATOM(J7))-1)+1)-LX(1)
                       YSHIFT=XYZ(NOPT*J4+3*(RBGROUP(NRUNNING-NRBSET+ALIGNATOM(J7))-1)+2)-LX(2)
                       ZSHIFT=XYZ(NOPT*J4+3*(RBGROUP(NRUNNING-NRBSET+ALIGNATOM(J7))-1)+3)-LX(3)
                    ELSE
                       XSHIFT=((NIMAGE-J4+1)*CMXS+J4*CMXF)/(NIMAGE+1)
                       YSHIFT=((NIMAGE-J4+1)*CMYS+J4*CMYF)/(NIMAGE+1)
                       ZSHIFT=((NIMAGE-J4+1)*CMZS+J4*CMZF)/(NIMAGE+1)
                    ENDIF
                    DO J5=1,NRBSET
                       IF (J4.LE.NIMAGE/2) THEN
                          LX(1:3)=LSTART(3*(J5-1)+1:3*(J5-1)+3)
                       ELSE
                          LX(1:3)=LFINISH(3*(J5-1)+1:3*(J5-1)+3)
                       ENDIF
                       CALL ROTATEANGLAXIS(LX,THETA,LVBEST)
                       LTEMP(3*(J5-1)+1)=LX(1)+XSHIFT
                       LTEMP(3*(J5-1)+2)=LX(2)+YSHIFT
                       LTEMP(3*(J5-1)+3)=LX(3)+ZSHIFT
                    ENDDO
                    DO J5=1,NRBSET
                       J6=RBGROUP(NRUNNING-NRBSET+J5)
                       COORDS(3*(J6-1)+1:3*(J6-1)+3)=LTEMP(3*(J5-1)+1:3*(J5-1)+3)
                    ENDDO
                    CALL POTENTIAL(COORDS,E1,LGDUMMY,.FALSE.,.FALSE.,RMSDUMMY,.FALSE.,.FALSE.)
!                   IF (DEBUG) PRINT '(A,G20.10,A,I6,A,G20.10)',' MAKEIMAGE> FOR ROTATION ',THETA,' IMAGE ',J4,' ENERGY=',E1
                    LEIMAGE2(J4)=E1
                    IF (E1.GT.EWORST2) THEN
                       JWORST2=J4
                       EWORST2=E1
                    ENDIF
                    IF (EWORST2.GT.EWORST*1.01D0) EXIT ! NO POINT CONTINUING - HIGHEST IMAGE IS HIGHER THAN PREVIOUS. 
                 ENDDO
                 IF (DEBUG) PRINT '(A,I6,A,G20.10)',' NNUTILS> HIGHEST IMAGE TESTED FOR + ROTATION WAS ',JWORST2, &
  &                                                 ' ENERGY=',EWORST2
!
!  TRY ROTATING THE OTHER WAY.
!
                 EWORST3=-HUGE(1.0D0)
                 IF (TBEST.LT.0.0D0) THEN
                    DUMMY=(6.283185307D0+TBEST)
                 ELSE
                    DUMMY=-(6.283185307D0-TBEST)
                 ENDIF
                 DO J4=1,NIMAGE
                    COORDS(1:NOPT)=XYZ(NOPT*J4+1:NOPT*J4+NOPT)

                    IF (J4.LE.NIMAGE/2) THEN
                       THETA=J4*DUMMY/(NIMAGE+1)
                    ELSE
                       THETA=-(NIMAGE-J4+1)*DUMMY/(NIMAGE+1)
                    ENDIF

                    IF (ALIGNCOMMON) THEN
                       IF (J4.LE.NIMAGE/2) THEN
                          LX(1:3)=LSTART(3*(ALIGNATOM(J7)-1)+1:3*(ALIGNATOM(J7)-1)+3)
                       ELSE
                          LX(1:3)=LFINISH(3*(ALIGNATOM(J7)-1)+1:3*(ALIGNATOM(J7)-1)+3)
                       ENDIF
                       CALL ROTATEANGLAXIS(LX,THETA,LVBEST)
                       XSHIFT=XYZ(NOPT*J4+3*(RBGROUP(NRUNNING-NRBSET+ALIGNATOM(J7))-1)+1)-LX(1)
                       YSHIFT=XYZ(NOPT*J4+3*(RBGROUP(NRUNNING-NRBSET+ALIGNATOM(J7))-1)+2)-LX(2)
                       ZSHIFT=XYZ(NOPT*J4+3*(RBGROUP(NRUNNING-NRBSET+ALIGNATOM(J7))-1)+3)-LX(3)
                    ELSE
                       XSHIFT=((NIMAGE-J4+1)*CMXS+J4*CMXF)/(NIMAGE+1)
                       YSHIFT=((NIMAGE-J4+1)*CMYS+J4*CMYF)/(NIMAGE+1)
                       ZSHIFT=((NIMAGE-J4+1)*CMZS+J4*CMZF)/(NIMAGE+1)
                    ENDIF
                    DO J5=1,NRBSET
                       IF (J4.LE.NIMAGE/2) THEN
                          LX(1:3)=LSTART(3*(J5-1)+1:3*(J5-1)+3)
                       ELSE
                          LX(1:3)=LFINISH(3*(J5-1)+1:3*(J5-1)+3)
                       ENDIF
                       CALL ROTATEANGLAXIS(LX,THETA,LVBEST)
                       LTEMP(3*(J5-1)+1)=LX(1)+XSHIFT
                       LTEMP(3*(J5-1)+2)=LX(2)+YSHIFT
                       LTEMP(3*(J5-1)+3)=LX(3)+ZSHIFT
                    ENDDO
                    DO J5=1,NRBSET
                       J6=RBGROUP(NRUNNING-NRBSET+J5)
                       COORDS(3*(J6-1)+1:3*(J6-1)+3)=LTEMP(3*(J5-1)+1:3*(J5-1)+3)
                    ENDDO
                    CALL POTENTIAL(COORDS,E2,LGDUMMY,.FALSE.,.FALSE.,RMSDUMMY,.FALSE.,.FALSE.)
!                   IF (DEBUG) PRINT '(A,G20.10,A,I6,A,G20.10)',' MAKEIMAGE> FOR ROTATION ',THETA,' IMAGE ',J4,' ENERGY=',E2
                    LEIMAGE3(J4)=E2
                    IF (E2.GT.MIN(EWORST3,EWORST*1.01D0)) THEN
                       JWORST3=J4
                       EWORST3=E2
                    ENDIF
                    IF (EWORST3.GT.MIN(EWORST,EWORST2)) EXIT ! NO POINT CONTINUING - HIGHEST IMAGE IS HIGHER THAN PREVIOUS. 
                 ENDDO
                 IF (DEBUG) PRINT '(A,I6,A,G20.10)',' NNUTILS> HIGHEST IMAGE TESTED FOR - ROTATION WAS ',JWORST3, &
  &                                                 ' ENERGY=',EWORST3

!
!  ALLOW THE ENERGY TO RISE BY ONE PER CENT?
!
                 IF (MIN(EWORST2,EWORST3).GT.EWORST*1.01D0) THEN
                       IF (DEBUG) PRINT '(A,G20.10,A)',' NNUTILS> BOTH INTERPOLATIONS RAISE THE ENERGY FROM ',LEIMAGE(JWORST),&
  &                                   ' - SKIP'
                       CYCLE ALIGNLOOP
                 ENDIF
                 IF (EWORST2.LT.EWORST3) THEN
                    LEIMAGE(1:NIMAGE)=LEIMAGE2(1:NIMAGE)
                    JWORST=JWORST2
                    DUMMY=TBEST
                    EWORST=EWORST2
                    IF (DEBUG) PRINT '(A,I6,A,I6,A)',' MAKEIMAGE> + ROTATIONAL INTERPOLATION FOR RIGID BODY GROUP ',J1, &
  &                                                  ' WITH ',NRBSET,' ATOMS'
                 ELSE
                    LEIMAGE(1:NIMAGE)=LEIMAGE3(1:NIMAGE)
                    JWORST=JWORST3
                    IF (TBEST.LT.0.0D0) THEN
                       DUMMY=(6.283185307D0+TBEST)
                    ELSE
                       DUMMY=-(6.283185307D0-TBEST)
                    ENDIF
                    EWORST=EWORST3
                    IF (DEBUG) PRINT '(A,I6,A,I6,A)',' MAKEIMAGE> - ROTATIONAL INTERPOLATION FOR RIGID BODY GROUP ',J1, &
  &                                                  ' WITH ',NRBSET,' ATOMS'
                 ENDIF

!                PRINT '(A,4F15.5)',' ANGLE, VECTOR=',DUMMY,LVBEST(1:3)
!                PRINT '(A)',' ATOM LIST:' 
!                PRINT '(22I6)',RBGROUP(NRUNNING-NRBSET+1:NRUNNING) 
!
!  NOW CHANGE THE POSITIONS OF THE ATOMS IN THIS GROUP IN EACH IMAGE, ROTATING THROUGH
!  ANGLE DUMMY EACH TIME. ROTATE LSTART FOR THE FIRST HALF, THEN LFINISH FOR THE 
!  SECOND HALF. SAME AS ABOVE, BUT WE ALREADY HAVE THE ENERGIES.
!
                 DO J4=1,NIMAGE
                    IF (J4.LE.NIMAGE/2) THEN
                       THETA=J4*DUMMY/(NIMAGE+1)
                    ELSE
                       THETA=-(NIMAGE-J4+1)*DUMMY/(NIMAGE+1)
                    ENDIF
                    IF (ALIGNCOMMON) THEN
                       IF (J4.LE.NIMAGE/2) THEN
                          LX(1:3)=LSTART(3*(ALIGNATOM(J7)-1)+1:3*(ALIGNATOM(J7)-1)+3)
                       ELSE
                          LX(1:3)=LFINISH(3*(ALIGNATOM(J7)-1)+1:3*(ALIGNATOM(J7)-1)+3)
                       ENDIF
                       CALL ROTATEANGLAXIS(LX,THETA,LVBEST)
                       XSHIFT=XYZ(NOPT*J4+3*(RBGROUP(NRUNNING-NRBSET+ALIGNATOM(J7))-1)+1)-LX(1)
                       YSHIFT=XYZ(NOPT*J4+3*(RBGROUP(NRUNNING-NRBSET+ALIGNATOM(J7))-1)+2)-LX(2)
                       ZSHIFT=XYZ(NOPT*J4+3*(RBGROUP(NRUNNING-NRBSET+ALIGNATOM(J7))-1)+3)-LX(3)
                    ELSE
                       XSHIFT=((NIMAGE-J4+1)*CMXS+J4*CMXF)/(NIMAGE+1)
                       YSHIFT=((NIMAGE-J4+1)*CMYS+J4*CMYF)/(NIMAGE+1)
                       ZSHIFT=((NIMAGE-J4+1)*CMZS+J4*CMZF)/(NIMAGE+1)
                    ENDIF
                    DO J5=1,NRBSET
                       IF (J4.LE.NIMAGE/2) THEN
                          LX(1:3)=LSTART(3*(J5-1)+1:3*(J5-1)+3)
                       ELSE
                          LX(1:3)=LFINISH(3*(J5-1)+1:3*(J5-1)+3)
                       ENDIF
                       CALL ROTATEANGLAXIS(LX,THETA,LVBEST)
                       LTEMP(3*(J5-1)+1)=LX(1)+XSHIFT
                       LTEMP(3*(J5-1)+2)=LX(2)+YSHIFT
                       LTEMP(3*(J5-1)+3)=LX(3)+ZSHIFT
                    ENDDO
                    DO J5=1,NRBSET
                       J6=RBGROUP(NRUNNING-NRBSET+J5)
                       ALIGNEDBEFORE(J6)=.TRUE.
                       XYZ(NOPT*J4+3*(J6-1)+1:NOPT*J4+3*(J6-1)+3)=LTEMP(3*(J5-1)+1:3*(J5-1)+3)
                    ENDDO
!!!! DJW
!                   COORDS(1:NOPT)=XYZ(NOPT*J4+1:NOPT*J4+NOPT)
!                   CALL POTENTIAL(COORDS,E2,LGDUMMY,.FALSE.,.FALSE.,RMSDUMMY,.FALSE.,.FALSE.)
!                   IF (DEBUG) PRINT '(A,I6,A,G20.10)',' MAKEIMAGE> CHECK IMAGE ',J4,' ENERGY=',E2
!!!! DJW
                 ENDDO
!
! DO NOT FORGET TO RESET THE END POINTS - THE RIGID BODY PERMUTATION MIGHT BE DIFFERENT!
!
                 DO J5=1,NRBSET
                    J6=RBGROUP(NRUNNING-NRBSET+J5)
                    XYZ(3*(J6-1)+1:3*(J6-1)+3)=QS(3*(J6-1)+1:3*(J6-1)+3)
                    XYZ((NIMAGE+1)*NOPT+3*(J6-1)+1:(NIMAGE+1)*NOPT+3*(J6-1)+3)=QF(3*(J6-1)+1:3*(J6-1)+3)
                 ENDDO
              ENDDO ALIGNLOOP
           ENDDO GROUPLOOP
           IF (NTRIES.LT.NRBTRIES) GOTO 11 ! DO IT AGAIN! DJW
        ENDIF
     ENDIF

     END SUBROUTINE MAKEIMAGE

     SUBROUTINE PRINTSUMMARY
          USE NEBDATA
          USE KEYNEB,ONLY: MOREPRINTING,NITERMAX
          USE CHARUTILS
          IMPLICIT NONE
          DOUBLE PRECISION :: TOTAL
          
          TOTAL = ENDTIME-STARTTIME
          IF (MOREPRINTING) THEN
               SELECT CASE (EXITSTATUS)
               CASE (1)
                    WRITE(*,'(1X,A)') 'CONVERGENCE CRITERION WAS SATISFIED.'
               CASE (2)
                    WRITE(*,'(1X,A)') 'REACHED MAXIMAL NUMBER OF ITERATIONS LIMIT.'
               END SELECT
               REALSTR=RM0S(WR(TOTAL,2))
               IF (NITERDONE.NE.0) THEN
                  REALSTR2=RM0S(WR(TOTAL/NITERDONE,2))
               ELSE
                  REALSTR2=RM0S(WR(TOTAL,2))
               ENDIF
               WRITE(*,'(1X,A)') 'TIME= '//TRIM(REALSTR)//' SEC ('//TRIM(REALSTR2)//' SEC/ITERATION)'
          ELSE ! THIS IS SHORT INFO WHICH IS A PART OF CONNECT OUTPUT
               IF (NITERMAX.GT.0) THEN
                  INTSTR = WI(NITERDONE)
                  WRITE(*,'(A)',ADVANCE='NO') ' DOUBLE-ENDED SEARCH ITERATIONS= '//TRIM(INTSTR)
                  REALSTR=WR(RMS,4)
                  WRITE(*,'(A)',ADVANCE='NO') ' RMS= '//TRIM(REALSTR)//' DEV= '
                  REALSTR=WR(AVDEV,2)
                  WRITE(*,'(A)',ADVANCE='NO') TRIM(REALSTR)//'% S= '
                  REALSTR=WR(SEPARATION,2)
                  WRITE(*,'(A)',ADVANCE='NO') TRIM(REALSTR)
                  REALSTR=WR(TOTAL,2)
                  WRITE(*,'(A)') ' TIME= '//TRIM(REALSTR)
               ENDIF
          ENDIF
     END SUBROUTINE PRINTSUMMARY

     SUBROUTINE DUMPFILES(I)
          USE NEBDATA
          USE PORFUNCS
          USE KEYNEB,ONLY: NIMAGE
          IMPLICIT NONE
          CHARACTER(LEN=1),INTENT(IN) :: I
          INTEGER ISTAT

          IF (I=="B") THEN
               OPEN(UNIT=90,FILE='RMSOFI',STATUS='REPLACE')
               OPEN(UNIT=91,FILE='AVDEVOFI',STATUS='REPLACE')
               OPEN(UNIT=92,FILE='EOFI',STATUS='REPLACE')
               OPEN(UNIT=93,FILE='SOFI',STATUS='REPLACE')
          ELSEIF (I=="E") THEN
               CLOSE(UNIT=90)
               CLOSE(UNIT=91)
               CLOSE(UNIT=92)
               CLOSE(UNIT=93)
          ELSEIF (I=="M") THEN
               WRITE(UNIT=90,FMT='(1X,I7,F20.10)') NITERDONE+NITERDONESAVE,RMS
               WRITE(UNIT=91,FMT='(1X,I7,F20.10)') NITERDONE+NITERDONESAVE,AVDEV
               WRITE(UNIT=92,FMT='(1X,I7,F20.10)') NITERDONE+NITERDONESAVE,ETOTAL/NIMAGE
               WRITE(UNIT=93,FMT='(1X,I7,F20.10)') NITERDONE+NITERDONESAVE,SEPARATION
               CALL FLUSH(90,ISTAT)
               CALL FLUSH(91,ISTAT)
               CALL FLUSH(92,ISTAT)
               CALL FLUSH(93,ISTAT)
          ENDIF
     END SUBROUTINE DUMPFILES

!    SUBROUTINE WRITEPROFILE(UNITIN)
     SUBROUTINE WRITEPROFILE(NITER)
          USE NEBDATA
          USE KEYNEB,ONLY: NIMAGE
          USE KEY,ONLY: FILTH,FILTHSTR
          IMPLICIT NONE

!         INTEGER,INTENT(IN),OPTIONAL :: UNITIN
          INTEGER,INTENT(IN) :: NITER
          
          INTEGER :: I,UNIT
          DOUBLE PRECISION :: DUMMY
          CHARACTER(LEN=20) :: FILENAME
          
!         IF (PRESENT(UNITIN)) THEN
!              UNIT=UNITIN
!         ELSE
               UNIT=992
               IF (NITER.GT.0) THEN
                  WRITE(FILENAME,'(I8)') NITER
                  FILENAME='NEB.EOFS.' // TRIM(ADJUSTL(FILENAME))
               ELSE   
                  FILENAME='NEB.EOFS'
               ENDIF
               IF (.NOT.FILTH==0) THEN
                    FILENAME=TRIM(FILENAME)//'.'//TRIM(ADJUSTL(FILTHSTR))
               ENDIF
               OPEN(UNIT=UNIT,FILE=FILENAME,STATUS='REPLACE')
!         ENDIF

          DUMMY=0.0D0
          WRITE(UNIT=UNIT,FMT='(2G24.13)') DUMMY,EEE(1)
          DO I=2,NIMAGE+1
               DUMMY = DUMMY + DVEC(I-1)
               WRITE(UNIT=UNIT,FMT='(2G24.13)') DUMMY,EEE(I)
          ENDDO
          DUMMY = DUMMY + DVEC(NIMAGE+1)
          WRITE(UNIT=UNIT,FMT='(2G24.13)') DUMMY,EEE(NIMAGE+2)

!         IF (.NOT.PRESENT(UNITIN)) THEN
               CLOSE(UNIT)
!         ENDIF
          PRINT *, 'WRITEPROFILE> NEB PROFILE WAS SAVED TO FILE "'//TRIM(FILENAME)//'"'
     END SUBROUTINE WRITEPROFILE

     SUBROUTINE RWG(WHAT,GUESS,NITER)
          USE PORFUNCS
          USE KEY,ONLY: FILTH,FILTHSTR,UNRST,STOCKT,AMHT,SEQ,NUMGLY,STOCKAAT, RBAAT,NTSITES
          USE COMMONS, ONLY: ZSYM, NRBSITES 
          USE INTCOMMONS, ONLY : DESMINT
          USE NEBDATA
          USE AMHGLOBALS, ONLY : NMRES
          USE KEYNEB,ONLY: NIMAGE,XYZFILE,RBXYZFILE,GUESSFILE
          IMPLICIT NONE

          CHARACTER,INTENT(IN) :: WHAT
          LOGICAL,INTENT(IN) :: GUESS
          INTEGER,INTENT(IN) :: NITER

          INTEGER :: EOF,J1,J2,J3,J4,GLY_COUNT
          CHARACTER(LEN=80) :: FILENAME,FILENAME2,DUMMYS,DUMMYS2
          DOUBLE PRECISION :: PEPCOORDS(3*NATOMS-6), STXYZ(3*NTSITES)

          IF (FILTH.EQ.0) THEN
             FILENAME=XYZFILE
             FILENAME2='NEB.PATH.UNR.XYZ'
             IF (GUESS) FILENAME=GUESSFILE
             IF (RBAAT) FILENAME2=RBXYZFILE
          ELSE
             FILENAME=TRIM(XYZFILE)//'.'//TRIM(ADJUSTL(FILTHSTR))
             FILENAME2='NEB.PATH.UNR.XYZ.'//TRIM(ADJUSTL(FILTHSTR))
             IF (GUESS) FILENAME=TRIM(GUESSFILE) !  //'.'//TRIM(ADJUSTL(FILTHSTR))
             IF (RBAAT) FILENAME2=TRIM(RBXYZFILE)//'.'//TRIM(ADJUSTL(FILTHSTR))
          ENDIF

          SELECT CASE(WHAT)
          CASE("W")
               IF (NITER.GT.0) THEN
                  IF (FILTH.EQ.0) THEN
                     WRITE(DUMMYS,'(I8)') NITER
                     DUMMYS2=TRIM(ADJUSTL(FILENAME))
                     FILENAME='NEB.' // TRIM(ADJUSTL(DUMMYS)) // '.XYZ' ! SO THAT VMD RECOGNISES THE FILE TYPE!
                     FILENAME2='RBNEB.' // TRIM(ADJUSTL(DUMMYS)) // '.XYZ'
                  ENDIF
               ENDIF
               OPEN(UNIT=993,FILE=FILENAME,STATUS='REPLACE')
               IF (DESMINT) THEN
                  DO J2=1,NIMAGE+2
                     WRITE(993,'(I4/)') NATOMS
                     WRITE(993,'(A5,1X,3F20.10)') (ZSYM((J1+2)/3),XYZCART( (J2-1)*NOPT+J1),&
                          & XYZCART((J2-1)*NOPT+J1+1), XYZCART((J2-1)*NOPT+J1+2),J1=1,NOPT,3)
                  ENDDO
               ELSEIF (STOCKT .OR. STOCKAAT) THEN
                  DO J2=1,NIMAGE+2
                     WRITE(993,'(I4/)') (NATOMS/2)
                     DO J1=1,(NATOMS/2)
                        WRITE(993,'(A5,1X,6F20.10)') ZSYM((J1+2)/3), &
                             & XYZ((J2-1)*NOPT+3*(J1-1)+1), XYZ((J2-1)*NOPT+3*(J1-1)+2), XYZ((J2-1)*NOPT+3*(J1-1)+3), &
  &    XYZ((J2-1)*NOPT+3*((NATOMS/2)+J1-1)+1), XYZ((J2-1)*NOPT+3*((NATOMS/2)+J1-1)+2), XYZ((J2-1)*NOPT+3*((NATOMS/2)+J1-1)+3)
                     ENDDO
                  ENDDO
               ELSEIF (RBAAT .AND. (.NOT. STOCKAAT)) THEN
                  OPEN(UNIT=114,FILE=FILENAME2,STATUS='UNKNOWN')
                  DO J2=1,NIMAGE+2
                     WRITE(993,'(I4/)') NATOMS/2
                     DO J1=1,(NATOMS/2)
                        WRITE(993,'(A5,1X,3F20.10)') 'O', &
                             & XYZ((J2-1)*NOPT+3*(J1-1)+1), XYZ((J2-1)*NOPT+3*(J1-1)+2), XYZ((J2-1)*NOPT+3*(J1-1)+3)
                     ENDDO
                     CALL SITEPOS(XYZ((J2-1)*NOPT+1:J2*NOPT),STXYZ)
                     WRITE(114,'(I4/)') (NATOMS/2)*NRBSITES
                     DO J1=1,(NATOMS/2)*NRBSITES
                        J3 = 3*J1
                        WRITE(114,'(A5,1X,3F20.10)') 'O', STXYZ(J3-2), STXYZ(J3-1), STXYZ(J3)
                     ENDDO
                  ENDDO
                  CLOSE(UNIT=114)
               ELSEIF (AMHT) THEN
                 DO J2=1,NIMAGE+2
!  GLY SET GETPARAMS.F 
!               WRITE(993,'(I4)')NATOMS +NUMGLY
!  GLY PRINTING TURNED OFF DJW
                WRITE(993,'(I4)')NATOMS 
                WRITE(993,*)'ENERGY'
                GLY_COUNT = 0
                
            DO J1=1,NMRES
              IF (SEQ(J1).EQ.8) THEN
          WRITE(993,'(A5,1X,3F20.10)') 'C1   ',XYZ((J2-1)*NOPT+9*(J1-1)+1-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+2-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+3-GLY_COUNT*3)
!  GLY PRINTING TURNED OFF DJW
!         WRITE(993,'(A5,1X,3F20.10)') 'C1   ',XYZ((J2-1)*NOPT+9*(J1-1)+1-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+2-GLY_COUNT*3), &
!    &                                  XYZ((J2-1)*NOPT+9*(J1-1)+3-GLY_COUNT*3)
          WRITE(993,'(A5,1X,3F20.10)') 'O    ',XYZ((J2-1)*NOPT+9*(J1-1)+4-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+5-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+6-GLY_COUNT*3)
                GLY_COUNT = GLY_COUNT +1
              ELSE
          WRITE(993,'(A5,1X,3F20.10)') 'C1   ',XYZ((J2-1)*NOPT+9*(J1-1)+1-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+2-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+3-GLY_COUNT*3)
          WRITE(993,'(A5,1X,3F20.10)') 'C2   ',XYZ((J2-1)*NOPT+9*(J1-1)+4-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+5-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+6-GLY_COUNT*3)
          WRITE(993,'(A5,1X,3F20.10)') 'O    ',XYZ((J2-1)*NOPT+9*(J1-1)+7-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+8-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+9-GLY_COUNT*3)
              ENDIF
          ENDDO

               ENDDO

                ELSE
                  DO J2=1,NIMAGE+2
                     WRITE(993,'(I4/)') NATOMS
                     WRITE(993,'(A5,1X,3F20.10)') (ZSYM((J1+2)/3),XYZ( (J2-1)*NOPT+J1),&
                          & XYZ((J2-1)*NOPT+J1+1), XYZ((J2-1)*NOPT+J1+2),J1=1,NOPT,3)
                  ENDDO
               ENDIF

               PRINT *, 'RWG> NEB COORDINATES WERE SAVED TO XYZ FILE "'//TRIM(FILENAME)//'"'

               IF (UNRST) THEN
                  OPEN(UNIT=114,FILE=FILENAME2,STATUS='UNKNOWN')
                  DO J2=1,NIMAGE+2
                       DO J3=1,(NATOMS/2)-1
                          DO J4=1,3
                             PEPCOORDS(6*(J3-1)+J4)=(2.0D0*XYZ((J2-1)*NOPT+(6*(J3-1)+J4))&
                             &+XYZ((J2-1)*NOPT+(6*J3+J4)))/3.0D0
                             PEPCOORDS(6*(J3-1)+J4+3)=(XYZ((J2-1)*NOPT+(6*(J3-1)+J4))+&
                             &2.0D0*XYZ((J2-1)*NOPT+(6*J3+J4)))/3.0D0
                          ENDDO
                       ENDDO
                       WRITE(114,'(I4/)') 2*NATOMS-2
                       WRITE(114,'(A5,1X,3F20.10)')('C    ',XYZ((J2-1)*NOPT+J1),&
                       &XYZ((J2-1)*NOPT+J1+1),XYZ((J2-1)*NOPT+J1+2),J1=1,NOPT,3)
                       WRITE(114,'(A5,1X,3F20.10)')('O    ', PEPCOORDS(J1),&
                       &PEPCOORDS(J1+1), PEPCOORDS(J1+2), J1=1,NOPT-11,6)
                       WRITE(114,'(A5,1X,3F20.10)')('N    ', PEPCOORDS(J1+3),&
                       &PEPCOORDS(J1+4), PEPCOORDS(J1+5), J1=1,NOPT-11,6)
                  ENDDO
                  CLOSE(UNIT=114)
                  PRINT *, 'COORDINATES WERE SAVED TO FILE ',FILENAME2
               ENDIF          
          CASE("R")
!              OPEN(UNIT=993,FILE=FILENAME,STATUS='OLD',IOSTAT=EOF)
               OPEN(UNIT=993,FILE=FILENAME,STATUS='OLD')
!              IF (.NOT.EOF==0) THEN
!                   PRINT *, "WHERE IS "//TRIM(ADJUSTL(FILENAME))//" FILE? - CAN'T FIND IT!"
!                   CALL TSUMMARY
!                   STOP
!              ENDIF
               DO J2=1,NIMAGE+2
!
! HERE WE SKIP TWO LINES, ALLOWING FOR THE SECOND LINE TO BE BLANK.
!
                     READ(993,'(A/)') DUMMYS 
                     READ(993,*) (ZSYM((J1+2)/3),XYZ( (J2-1)*NOPT+J1),&
                          & XYZ((J2-1)*NOPT+J1+1), XYZ((J2-1)*NOPT+J1+2),J1=1,NOPT,3)
               ENDDO
               PRINT '(A)','NNUTILS> GUESS READ FROM FILE ' // TRIM(ADJUSTL(GUESSFILE))
          END SELECT
          CLOSE(UNIT=993)
     END SUBROUTINE RWG

     SUBROUTINE SAVEBANDCOORD
          USE NEBDATA,ONLY: XYZ,NOPT
          USE KEYNEB,ONLY:NIMAGE,PTSFILE
          IMPLICIT NONE

          INTEGER :: RECLEN,I

          INQUIRE(IOLENGTH=RECLEN) XYZ(1:NOPT)
          OPEN(UNIT=40,FILE=TRIM(PTSFILE),STATUS='UNKNOWN',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=RECLEN)
          DO I=1,NIMAGE+2
               WRITE(40,REC=I) ( XYZ(NOPT*(I-1)+1:NOPT*I) )
          ENDDO
          CLOSE(40)
          PRINT *, 'NEB COORDINATES WERE SAVED TO BINARY FILE "'//TRIM(PTSFILE)//'"'
     END SUBROUTINE SAVEBANDCOORD
END MODULE NEBUTILS
