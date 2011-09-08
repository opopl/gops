!   NEB module is an implementation of the nudged elastic band method for performing double-ended pathway searches.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of NEB module. NEB module is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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

          ! /------------------Separation(along the MEP)-------------------\
          ! Q--------Image1------Image2---- ...----ImageNimage-------------Fin
          ! \_DVec(1)_/ \_DVec(2)_/                       \_DVec(Nimage+1)_/
     
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
             PRINT '(A)','internalimagedistribution> ERROR - this routine should not be called with BULKT true'
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
           
          ! for internals
          DOUBLE PRECISION :: DELTACART(3*NATOMS)
          LOGICAL :: INTPTEST, FAILED
          ! alignment stuff
          DOUBLE PRECISION :: DISTF, DIST, DIST2, RMAT(3,3)
          CHARACTER(LEN=5) :: ZSYMSAVE
          COMMON /SYS/ ZSYMSAVE

          IF (DESMINT) THEN
             BACKTCUTOFF = INTERPBACKTCUT
             INTPTEST = .FALSE.
          ENDIF

          IF (MORPHT) THEN
             IF (DESMINT) THEN
                print*, 'ERROR! MORPHT not implemented with DESM.'
                STOP
             ENDIF

             QS(1:3*NATOMS)=XYZ(1:NOPT)
             QF(1:3*NATOMS)=XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+1)+NOPT)
             KNOWE=.FALSE.;  KNOWG=.FALSE.; KNOWH=.FALSE.
          
             CALL MORPH(MSTEPS,QS,QF,ENERGY,VNEW,MFLAG,LRMS,ITDONE,.TRUE.)
             IF (ITDONE.LT.NIMAGE) THEN
                IF (MFLAG) THEN
!
!  Unfortunately many statements in the NEB
!  routines do not specify the array bounds on array operations, which then default to the
!  declared array size and therefore go wrong. Hence if we reduce the number of images here
!  then we have to deallocate and reallocate:
!                  
                   PRINT '(A,I6)',' makeimage> Fewer morph points than images - reducing images to ',ITDONE
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

                   INTERVAL=ITDONE/Nimage
                   NDONE=0
                   OPEN(UNIT=991,FILE='morph.points',STATUS='OLD')
                   DO J1=2,NIMAGE+1
                      DO J2=1,INTERVAL
                         NDONE=NDONE+1
                         READ(991,*) xyz(nopt*(J1-1)+1:nopt*(J1-1)+NOPT) 
                      ENDDO
                      IF (DEBUG) PRINT '(2(A,I6))',' makeimage> Image ',J1,' read from morph.points frame ',NDONE
                   ENDDO
                   CLOSE(991)
                   xyz(nopt*(Nimage+1)+1:nopt*(Nimage+1)+NOPT)=QF(1:NOPT) ! to align the final image
                ELSE
                   PRINT '(A)',' makeimage> Fewer morph points than images - revert to linear interpolation'
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
                OPEN(UNIT=991,FILE='morph.points',STATUS='OLD')
                DO J1=2,NIMAGE+1
                   DO J2=1,INTERVAL
                      NDONE=NDONE+1
                      READ(991,*) XYZ(NOPT*(J1-1)+1:NOPT*(J1-1)+NOPT) 
                   ENDDO
                   IF (DEBUG) PRINT '(2(A,I6))',' makeimage> Image ',J1,' read from morph.points frame ',NDONE
                ENDDO
                CLOSE(991)
                XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+1)+NOPT)=QF(1:NOPT) ! TO ALIGN THE FINAL IMAGE
             ENDIF
             MORPHT=.FALSE. ! DJW
          ELSEIF (GREATCIRCLET) THEN 
             IF (NOPT.NE.3*NATOMS) THEN
                PRINT '(A)','nnutils> ERROR - NOPT needs to be 3*NATOMS to use GREATCIRCLE interpolation'
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

!     QUATERNION INTERPOLATION: iSLERP
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
   
! bs360: Interpolation using dihedrals.
             IF (CHRMMT.AND.CHICDNEB.AND.ICINTERPT) THEN
                QS(1:NOPT)=XYZ(1:NOPT)         
                QF(1:NOPT)=XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+2))
             ENDIF
! end bs360
! msb50: Interpolation using dihedrals for amber
             IF ((AMBERT.OR.NABT).AND.AMBICDNEBT.AND.AMBERICT) THEN
                QS(1:NOPT)=XYZ(1:NOPT)
                QF(1:NOPT)=XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+2))
                PRINT*, "msb50 NOPT", NOPT
!                DO K=1, 3*NATOMS
!                    PRINT '(6f9.3)', QF(6*(K-1)+1:6*K)
!                ENDDO
             ENDIF
! end msb50
             DO I=2,NIMAGE+1
                  XYZ(NOPT*(I-1)+1:NOPT*I) = XYZ(1:NOPT) + DELTAX*(I-1)
       
! bs360
                  IF (CHRMMT.AND.CHICDNEB.AND.ICINTERPT) THEN
                     COORDS(1:NOPT)=XYZ(NOPT*(I-1)+1:NOPT*I)
                     SFRAC=1.D0-(I-1.D0)/(NIMAGE+1.D0)                  
!                     write(6,*) 'i,NIMAGE+1,sfrac= ',i,NIMAGE+1,sfrac
                     CALL ICINTERPOL(COORDS,QS,QF,SFRAC)
                     XYZ(NOPT*(I-1)+1:NOPT*I)=COORDS(1:NOPT)        
!                     write(100+I,'(3F10.3)') XYZ(NOPT*(I-1)+1:NOPT*I)
!                     call flush(100+I)
                     call flush(100+I,ISTAT)
                  ENDIF
! end bs360

! msb50
                  IF ((AMBERT.OR.NABT).AND.AMBICDNEBT.AND.AMBERICT) THEN
                     COORDS(1:NOPT)=XYZ(NOPT*(I-1)+1:NOPT*I)
                     SFRAC=1.D0-(I-1.D0)/(NIMAGE+1.D0)
                     IF (.NOT.ALLOCATED(NICTOT)) THEN
                         CALL SETDIHEAM() 
                     ENDIF
                     CALL TAKESTEPAMDIHED(COORDS, QS, QF,SFRAC)
                     write(6,*) 'i,NIMAGE+1,sfrac= ',i,NIMAGE+1,sfrac
                     XYZ(NOPT*(I-1)+1:NOPT*I)=COORDS(1:NOPT)
                     !DO k=1,(NOPT/3)
                     !  PRINT*, XYZ(3*(k-1)+1:3*k) 
                     !ENDDO
!                    write(100+I,'(3F10.3)') XYZ(NOPT*(I-1)+1:NOPT*I)
!                     call flush(100+I)
                  ENDIF
! end msb50

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
!  Now if we have some rigid bodies, interpolate them as rigid bodies!
!  However, also check whether the alternative interpolation based on
!  regular minpermdist for the whole molecule gives a better initial guess.
!  If so, we can replace blocks of atoms with this interpolation, but
!  we must use the same local permutational isomer for each image.
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
!  Identify the current highest image. We will use this to choose between
!  + and - rotations about the chosen RB axis for each image.
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
!             IF (DEBUG) PRINT '(A,I6,A,G20.10)',' nnutils> image ',J1,' energy=',E1
           ENDDO
           IF (DEBUG) PRINT '(A,I6,A,G20.10)',' nnutils> highest image is ',JWORST,' energy=',EWORST
           NTRIES=NTRIES+1
!          PRINT '(A,I6)',' makeimage> number of entries in RBGROUP=',NRBTOTAL
!          PRINT '(22I6)',RBGROUP(1:NRBTOTAL)
!
!  Rigid body interpolations use LSTART and LFINISH, which are assigned from QS and QF.
!  We replace blocks of atoms in the XYZ array corresponding to these atoms if the local
!  rigid body interpolation lowers the energy.
!
           QS(1:3*NATOMS)=XYZ(1:NOPT)
           QF(1:3*NATOMS)=XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+1)+NOPT)
!
!  Now change the whole XYZ array to be based on the permutational isomer that
!  minimises the overall distance. The difference is in the initial straight line
!  guess, which may be for different permutational isomers.
!  Use the permutation that gives the lowest value for the highest image.
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
!             IF (DEBUG) PRINT '(A,I6,A,G20.10)',' nnutils> alternative interpolation, image ',J1,' energy=',E1
              IF (EWORST2.GT.EWORST) EXIT ! no point continuing, we will use the other initial guess.
           ENDDO
           IF (EWORST2.LT.EWORST) THEN
              IF (DEBUG) PRINT '(A)',' nnutils> adopting alternative interpolation for initial guess'
              JWORST=JWORST2
              EWORST=EWORST2 
              LEIMAGE(1:NIMAGE)=LEIMAGE2(1:NIMAGE)
              IF (DEBUG) PRINT '(A,I6,A,G20.10)',' nnutils> highest image is ',JWORST,' energy=',EWORST
              XYZ(1:3*NATOMS)=QS2(1:3*NATOMS)
              XYZ((NIMAGE+1)*NOPT+1:(NIMAGE+1)*NOPT+3*NATOMS)=QF2(1:3*NATOMS)
              DO J1=1,NIMAGE
                 XYZ(J1*NOPT+1:J1*NOPT+3*NATOMS)=(J1*QF2(1:3*NATOMS)+(NIMAGE+1-J1)*QS2(1:3*NATOMS))/(NIMAGE+1)
              ENDDO
           ENDIF
!
! Initial guess has been selected
!
           NRUNNING=0
           grouploop: DO J1=1,NRBGROUP
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
!                PRINT '(3(A,I6))',' makeimage> group ',J1,' contains atom ',J3,' total in group=',RBNINGROUP(J1)
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
                 PRINT '(A,I8,A,I8)',' makeimage> Number of atoms in rigid body group ',J1,' is ',NRBSET
!                PRINT '(A,3F15.5)',' makeimage> start  centre of mass: ',CMXS,CMYS,CMZS
!                PRINT '(A,3F15.5)',' makeimage> start  coordinates:'
!                PRINT '(3F20.10)',LSTART(1:3*NRBSET)
!                PRINT '(A,3F15.5)',' makeimage> finish centre of mass: ',CMXF,CMYF,CMZF
!                PRINT '(A,3F15.5)',' makeimage> finish coordinates:'
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
! Check result by back-rotation and find the best angle-axis that gives the minimum error.
!
! Dot product in DUMMY should be the same for the atom in start and finish.
! However, we are dealing with approximate reigid bodies, so the symmetry could be broken
! with a likely loss of accuracy.
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
! Try theta1
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
!                   IF (DEBUG) PRINT '(A,2F20.10,A,3F15.5)',' makeimage> theta1, residual: ', &
! &                            THETA1,SQRT(DUMMY),' axis: ',LV(1:3)  ! DJW
!
! Try theta2
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
!                   IF (DEBUG) PRINT '(A,2F20.10,A,3F15.5)',' makeimage> theta2, residual: ', &
! &                            THETA2,SQRT(DUMMY),' axis: ',LV(1:3)  ! DJW
                    THETA1=-THETA1; THETA2=-THETA2
!
! Try -theta1
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
!                   IF (DEBUG) PRINT '(A,2F20.10,A,3F15.5)',' makeimage> theta1, residual: ', &
! &                            THETA1,SQRT(DUMMY),' axis: ',LV(1:3)  ! DJW
!
! Try -theta2
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
!                   IF (DEBUG) PRINT '(A,2F20.10,A,3F15.5)',' makeimage> theta2, residual: ', &
! &                            THETA2,SQRT(DUMMY),' axis: ',LV(1:3)  ! DJW
                 ENDDO
              ENDDO
!
! If the angle is small (radians) then presumably a straight line interpolation should be OK.
! TBEST is also set to zero if the attempted rigid body rotation was not accurate enough.
!
              IF (ABS(TBEST).LT.0.1D0) CYCLE grouploop

              IF (DEBUG) THEN
                 PRINT '(2(A,G20.10))',' makeimage> best residual is ',SQRT(DBEST),' for angle ',TBEST
!                PRINT '(A,3F15.5)',' makeimage> rotated coordinates:'
!                DO J4=1,NRBSET
!                   LX(1:3)=LSTART(3*(J4-1)+1:3*(J4-1)+3)
!                   CALL ROTATEANGLAXIS(LX,TBEST,LVBEST)
!                   DUMMY=DUMMY+(LX(1)-LFINISH(3*(J4-1)+1))**2+(LX(2)-LFINISH(3*(J4-1)+2))**2+(LX(3)-LFINISH(3*(J4-1)+3))**2
!                   PRINT '(3F20.10)',LX(1:3)
!                ENDDO
              ENDIF

!
!  Interpolate this group using the best rotation. 
!  We can try TBEST/2 and (2*pi-TBEST)/2 to try and guess the best path.
!
!  If the group contains an atom that is in a previous rigid body then translate the
!  group so that this atom lines up with its previous images. Otherwise shift the
!  centre of mass sequentially between the values for the starting and finishing 
!  geometries. It is possible that the current group will contain more than one atom
!  in common with previous groups. In this case we use the last atom.
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
                       IF (DEBUG) PRINT '(A,I6,A)',' makeimage> atom ',RBGROUP(NRUNNING-NRBSET+J5), &
  &                                                ' was in a previously aligned rigid group'
                    ENDIF
                 ENDDO
              ENDDO
              IF (DEBUG.AND.ALIGNCOMMON) PRINT '(A,I6,A)',' makeimage> translating group to align each shared atom in turn '

              alignloop: DO J7=1,MAX(NALIGNATOM,1)
                 IF (DEBUG.AND.ALIGNCOMMON) PRINT '(A,I6,A)',' makeimage> aligning on atom ',RBGROUP(NRUNNING-NRBSET+ALIGNATOM(J7))

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
!                   IF (DEBUG) PRINT '(A,G20.10,A,I6,A,G20.10)',' makeimage> for rotation ',THETA,' image ',J4,' energy=',E1
                    LEIMAGE2(J4)=E1
                    IF (E1.GT.EWORST2) THEN
                       JWORST2=J4
                       EWORST2=E1
                    ENDIF
                    IF (EWORST2.GT.EWORST*1.01D0) EXIT ! no point continuing - highest image is higher than previous. 
                 ENDDO
                 IF (DEBUG) PRINT '(A,I6,A,G20.10)',' nnutils> highest image tested for + rotation was ',JWORST2, &
  &                                                 ' energy=',EWORST2
!
!  Try rotating the other way.
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
!                   IF (DEBUG) PRINT '(A,G20.10,A,I6,A,G20.10)',' makeimage> for rotation ',THETA,' image ',J4,' energy=',E2
                    LEIMAGE3(J4)=E2
                    IF (E2.GT.MIN(EWORST3,EWORST*1.01D0)) THEN
                       JWORST3=J4
                       EWORST3=E2
                    ENDIF
                    IF (EWORST3.GT.MIN(EWORST,EWORST2)) EXIT ! no point continuing - highest image is higher than previous. 
                 ENDDO
                 IF (DEBUG) PRINT '(A,I6,A,G20.10)',' nnutils> highest image tested for - rotation was ',JWORST3, &
  &                                                 ' energy=',EWORST3

!
!  Allow the energy to rise by one per cent?
!
                 IF (MIN(EWORST2,EWORST3).GT.EWORST*1.01D0) THEN
                       IF (DEBUG) PRINT '(A,G20.10,A)',' nnutils> both interpolations raise the energy from ',LEIMAGE(JWORST),&
  &                                   ' - skip'
                       CYCLE alignloop
                 ENDIF
                 IF (EWORST2.LT.EWORST3) THEN
                    LEIMAGE(1:NIMAGE)=LEIMAGE2(1:NIMAGE)
                    JWORST=JWORST2
                    DUMMY=TBEST
                    EWORST=EWORST2
                    IF (DEBUG) PRINT '(A,I6,A,I6,A)',' makeimage> + rotational interpolation for rigid body group ',J1, &
  &                                                  ' with ',NRBSET,' atoms'
                 ELSE
                    LEIMAGE(1:NIMAGE)=LEIMAGE3(1:NIMAGE)
                    JWORST=JWORST3
                    IF (TBEST.LT.0.0D0) THEN
                       DUMMY=(6.283185307D0+TBEST)
                    ELSE
                       DUMMY=-(6.283185307D0-TBEST)
                    ENDIF
                    EWORST=EWORST3
                    IF (DEBUG) PRINT '(A,I6,A,I6,A)',' makeimage> - rotational interpolation for rigid body group ',J1, &
  &                                                  ' with ',NRBSET,' atoms'
                 ENDIF

!                PRINT '(A,4F15.5)',' angle, vector=',DUMMY,LVBEST(1:3)
!                PRINT '(A)',' atom list:' 
!                PRINT '(22I6)',RBGROUP(NRUNNING-NRBSET+1:NRUNNING) 
!
!  Now change the positions of the atoms in this group in each image, rotating through
!  angle DUMMY each time. Rotate LSTART for the first half, then LFINISH for the 
!  second half. Same as above, but we already have the energies.
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
!                   IF (DEBUG) PRINT '(A,I6,A,G20.10)',' makeimage> CHECK image ',J4,' energy=',E2
!!!! DJW
                 ENDDO
!
! Do not forget to reset the end points - the rigid body permutation might be different!
!
                 DO J5=1,NRBSET
                    J6=RBGROUP(NRUNNING-NRBSET+J5)
                    XYZ(3*(J6-1)+1:3*(J6-1)+3)=QS(3*(J6-1)+1:3*(J6-1)+3)
                    XYZ((NIMAGE+1)*NOPT+3*(J6-1)+1:(NIMAGE+1)*NOPT+3*(J6-1)+3)=QF(3*(J6-1)+1:3*(J6-1)+3)
                 ENDDO
              ENDDO alignloop
           ENDDO grouploop
           IF (NTRIES.LT.NRBTRIES) GOTO 11 ! do it again! DJW
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
                    WRITE(*,'(1x,a)') 'Convergence criterion was satisfied.'
               CASE (2)
                    WRITE(*,'(1x,a)') 'Reached maximal number of iterations limit.'
               END SELECT
               REALSTR=RM0S(WR(TOTAL,2))
               IF (NITERDONE.NE.0) THEN
                  REALSTR2=RM0S(WR(TOTAL/NITERDONE,2))
               ELSE
                  REALSTR2=RM0S(WR(TOTAL,2))
               ENDIF
               WRITE(*,'(1x,a)') 'time= '//trim(RealStr)//' sec ('//trim(RealStr2)//' sec/iteration)'
          ELSE ! THIS IS SHORT INFO WHICH IS A PART OF CONNECT OUTPUT
               IF (NITERMAX.GT.0) THEN
                  INTSTR = WI(NITERDONE)
                  WRITE(*,'(a)',advance='no') ' Double-ended search iterations= '//trim(IntStr)
                  REALSTR=WR(RMS,4)
                  WRITE(*,'(a)',advance='no') ' RMS= '//trim(RealStr)//' Dev= '
                  REALSTR=WR(AVDEV,2)
                  WRITE(*,'(a)',advance='no') trim(RealStr)//'% S= '
                  REALSTR=WR(SEPARATION,2)
                  WRITE(*,'(a)',advance='no') trim(RealStr)
                  REALSTR=WR(TOTAL,2)
                  WRITE(*,'(a)') ' time= '//trim(RealStr)
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

          IF (I=="b") then
               OPEN(UNIT=90,FILE='RmsofI',status='replace')
               OPEN(UNIT=91,FILE='AvDevofI',status='replace')
               OPEN(UNIT=92,FILE='EofI',status='replace')
               OPEN(UNIT=93,FILE='SofI',status='replace')
          ELSEIF (I=="e") then
               CLOSE(UNIT=90)
               CLOSE(UNIT=91)
               CLOSE(UNIT=92)
               CLOSE(UNIT=93)
          ELSEIF (I=="m") then
               WRITE(UNIT=90,FMT='(1x,i7,f20.10)') NIterDone+NIterDoneSave,RMS
               WRITE(UNIT=91,FMT='(1x,i7,f20.10)') NIterDone+NIterDoneSave,AvDev
               WRITE(UNIT=92,FMT='(1x,i7,f20.10)') NIterDone+NIterDoneSave,Etotal/Nimage
               WRITE(UNIT=93,FMT='(1x,i7,f20.10)') NIterDone+NIterDoneSave,Separation
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
                  FILENAME='neb.EofS.' // TRIM(ADJUSTL(FILENAME))
               ELSE   
                  FILENAME='neb.EofS'
               ENDIF
               IF (.NOT.FILTH==0) THEN
                    FILENAME=TRIM(FILENAME)//'.'//TRIM(ADJUSTL(FILTHSTR))
               ENDIF
               OPEN(UNIT=UNIT,FILE=FILENAME,STATUS='replace')
!         ENDIF

          DUMMY=0.0D0
          WRITE(UNIT=UNIT,FMT='(2g24.13)') dummy,eee(1)
          DO I=2,NIMAGE+1
               DUMMY = DUMMY + DVEC(I-1)
               WRITE(UNIT=UNIT,FMT='(2g24.13)') dummy,eee(i)
          ENDDO
          DUMMY = DUMMY + DVEC(NIMAGE+1)
          WRITE(UNIT=UNIT,FMT='(2g24.13)') dummy,eee(Nimage+2)

!         IF (.NOT.PRESENT(UNITIN)) THEN
               CLOSE(UNIT)
!         ENDIF
          PRINT *, 'writeprofile> NEB profile was saved to file "'//trim(filename)//'"'
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
             FILENAME2='neb.path.unr.xyz'
             IF (GUESS) FILENAME=GUESSFILE
             IF (RBAAT) FILENAME2=RBXYZFILE
          ELSE
             FILENAME=TRIM(XYZFILE)//'.'//TRIM(ADJUSTL(FILTHSTR))
             FILENAME2='neb.path.unr.xyz.'//TRIM(ADJUSTL(FILTHSTR))
             IF (GUESS) FILENAME=TRIM(GUESSFILE) !  //'.'//TRIM(ADJUSTL(FILTHSTR))
             IF (RBAAT) FILENAME2=TRIM(RBXYZFILE)//'.'//TRIM(ADJUSTL(FILTHSTR))
          ENDIF

          SELECT CASE(WHAT)
          CASE("w")
               IF (NITER.GT.0) THEN
                  IF (FILTH.EQ.0) THEN
                     WRITE(DUMMYS,'(I8)') NITER
                     DUMMYS2=TRIM(ADJUSTL(FILENAME))
                     FILENAME='neb.' // TRIM(ADJUSTL(DUMMYS)) // '.xyz' ! so that vmd recognises the file type!
                     FILENAME2='rbneb.' // TRIM(ADJUSTL(DUMMYS)) // '.xyz'
                  ENDIF
               ENDIF
               OPEN(UNIT=993,FILE=FILENAME,STATUS='replace')
               IF (DESMINT) THEN
                  DO J2=1,NIMAGE+2
                     WRITE(993,'(i4/)') natoms
                     WRITE(993,'(a5,1x,3f20.10)') (ZSYM((j1+2)/3),XYZCART( (j2-1)*Nopt+j1),&
                          & XYZCART((J2-1)*NOPT+J1+1), XYZCART((J2-1)*NOPT+J1+2),J1=1,NOPT,3)
                  ENDDO
               ELSEIF (STOCKT .OR. STOCKAAT) THEN
                  DO J2=1,NIMAGE+2
                     WRITE(993,'(i4/)') (natoms/2)
                     DO J1=1,(natoms/2)
                        WRITE(993,'(a5,1x,6f20.10)') ZSYM((j1+2)/3), &
                             & XYZ((J2-1)*NOPT+3*(J1-1)+1), XYZ((J2-1)*NOPT+3*(J1-1)+2), XYZ((J2-1)*NOPT+3*(J1-1)+3), &
  &    XYZ((J2-1)*NOPT+3*((natoms/2)+J1-1)+1), XYZ((J2-1)*NOPT+3*((natoms/2)+J1-1)+2), XYZ((J2-1)*NOPT+3*((natoms/2)+J1-1)+3)
                     ENDDO
                  ENDDO
               ELSEIF (RBAAT .AND. (.NOT. STOCKAAT)) THEN
                  OPEN(UNIT=114,FILE=FILENAME2,STATUS='unknown')
                  DO J2=1,NIMAGE+2
                     WRITE(993,'(i4/)') NATOMS/2
                     DO J1=1,(NATOMS/2)
                        WRITE(993,'(a5,1x,3f20.10)') 'O', &
                             & XYZ((J2-1)*NOPT+3*(J1-1)+1), XYZ((J2-1)*NOPT+3*(J1-1)+2), XYZ((J2-1)*NOPT+3*(J1-1)+3)
                     ENDDO
                     CALL SITEPOS(XYZ((J2-1)*NOPT+1:J2*NOPT),STXYZ)
                     WRITE(114,'(i4/)') (NATOMS/2)*NRBSITES
                     DO J1=1,(NATOMS/2)*NRBSITES
                        J3 = 3*J1
                        WRITE(114,'(a5,1x,3f20.10)') 'O', STXYZ(J3-2), STXYZ(J3-1), STXYZ(J3)
                     ENDDO
                  ENDDO
                  CLOSE(UNIT=114)
               ELSEIF (AMHT) THEN
                 DO J2=1,NIMAGE+2
!  GLY set getparams.f 
!               WRITE(993,'(i4)')NATOMS +NUMGLY
!  GLY printing turned off DJW
                WRITE(993,'(i4)')NATOMS 
                WRITE(993,*)'Energy'
                GLY_COUNT = 0
                
            DO J1=1,NMRES
              IF (SEQ(J1).EQ.8) THEN
          WRITE(993,'(a5,1x,3f20.10)') 'C1   ',XYZ((J2-1)*NOPT+9*(J1-1)+1-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+2-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+3-GLY_COUNT*3)
!  GLY printing turned off DJW
!         WRITE(993,'(a5,1x,3f20.10)') 'C1   ',XYZ((J2-1)*NOPT+9*(J1-1)+1-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+2-GLY_COUNT*3), &
!    &                                  XYZ((J2-1)*NOPT+9*(J1-1)+3-GLY_COUNT*3)
          WRITE(993,'(a5,1x,3f20.10)') 'O    ',XYZ((J2-1)*NOPT+9*(J1-1)+4-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+5-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+6-GLY_COUNT*3)
                GLY_COUNT = GLY_COUNT +1
              ELSE
          WRITE(993,'(a5,1x,3f20.10)') 'C1   ',XYZ((J2-1)*NOPT+9*(J1-1)+1-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+2-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+3-GLY_COUNT*3)
          WRITE(993,'(a5,1x,3f20.10)') 'C2   ',XYZ((J2-1)*NOPT+9*(J1-1)+4-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+5-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+6-GLY_COUNT*3)
          WRITE(993,'(a5,1x,3f20.10)') 'O    ',XYZ((J2-1)*NOPT+9*(J1-1)+7-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+8-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+9-GLY_COUNT*3)
              ENDIF
          ENDDO

               ENDDO

                ELSE
                  DO J2=1,NIMAGE+2
                     WRITE(993,'(i4/)') natoms
                     WRITE(993,'(a5,1x,3f20.10)') (ZSYM((j1+2)/3),xyz( (j2-1)*Nopt+j1),&
                          & XYZ((J2-1)*NOPT+J1+1), XYZ((J2-1)*NOPT+J1+2),J1=1,NOPT,3)
                  ENDDO
               ENDIF

               PRINT *, 'rwg> NEB coordinates were saved to xyz file "'//trim(filename)//'"'

               IF (UNRST) THEN
                  OPEN(UNIT=114,FILE=FILENAME2,STATUS='unknown')
                  DO J2=1,NIMAGE+2
                       DO J3=1,(NATOMS/2)-1
                          DO J4=1,3
                             PEPCOORDS(6*(J3-1)+J4)=(2.0D0*XYZ((J2-1)*NOPT+(6*(J3-1)+J4))&
                             &+XYZ((J2-1)*NOPT+(6*J3+J4)))/3.0D0
                             PEPCOORDS(6*(J3-1)+J4+3)=(XYZ((J2-1)*NOPT+(6*(J3-1)+J4))+&
                             &2.0D0*XYZ((J2-1)*NOPT+(6*J3+J4)))/3.0D0
                          ENDDO
                       ENDDO
                       WRITE(114,'(i4/)') 2*natoms-2
                       WRITE(114,'(a5,1x,3f20.10)')('C    ',xyz((j2-1)*Nopt+j1),&
                       &XYZ((J2-1)*NOPT+J1+1),XYZ((J2-1)*NOPT+J1+2),J1=1,NOPT,3)
                       WRITE(114,'(a5,1x,3f20.10)')('O    ', pepcoords(j1),&
                       &PEPCOORDS(J1+1), PEPCOORDS(J1+2), J1=1,NOPT-11,6)
                       WRITE(114,'(a5,1x,3f20.10)')('N    ', pepcoords(j1+3),&
                       &PEPCOORDS(J1+4), PEPCOORDS(J1+5), J1=1,NOPT-11,6)
                  ENDDO
                  CLOSE(UNIT=114)
                  PRINT *, 'Coordinates were saved to file ',filename2
               ENDIF          
          CASE("r")
!              OPEN(UNIT=993,FILE=FILENAME,STATUS='old',iostat=eof)
               OPEN(UNIT=993,FILE=FILENAME,STATUS='old')
!              IF (.NOT.EOF==0) THEN
!                   PRINT *, "WHERE IS "//TRIM(ADJUSTL(FILENAME))//" FILE? - CAN't find it!"
!                   CALL TSUMMARY
!                   STOP
!              ENDIF
               DO J2=1,NIMAGE+2
!
! Here we skip two lines, allowing for the second line to be blank.
!
                     READ(993,'(A/)') DUMMYS 
                     READ(993,*) (ZSYM((j1+2)/3),XYZ( (J2-1)*NOPT+J1),&
                          & XYZ((J2-1)*NOPT+J1+1), XYZ((J2-1)*NOPT+J1+2),J1=1,NOPT,3)
               ENDDO
               PRINT '(A)','nnutils> Guess read from file ' // TRIM(ADJUSTL(GUESSFILE))
          END SELECT
          CLOSE(UNIT=993)
     END SUBROUTINE RWG

     SUBROUTINE SAVEBANDCOORD
          USE NEBDATA,ONLY: XYZ,NOPT
          USE KEYNEB,ONLY:NIMAGE,PTSFILE
          IMPLICIT NONE

          INTEGER :: RECLEN,I

          INQUIRE(IOLENGTH=RECLEN) XYZ(1:NOPT)
          OPEN(UNIT=40,FILE=TRIM(PTSFILE),STATUS='unknown',form='unformatted',access='direct',recl=reclen)
          DO I=1,NIMAGE+2
               WRITE(40,REC=I) ( XYZ(NOPT*(I-1)+1:NOPT*I) )
          ENDDO
          CLOSE(40)
          PRINT *, 'NEB coordinates were saved to binary file "'//trim(PtsFile)//'"'
     END SUBROUTINE SAVEBANDCOORD
END MODULE NEBUTILS
