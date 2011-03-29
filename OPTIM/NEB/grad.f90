!   NEB MODULE IS AN IMPLEMENTATI ON OF THE NUDGED ELASTIC BAND METHOD FOR PERFORMING DOUBLE-ENDED PATHWAY SEARCHES.
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
MODULE GRADIENTS
     IMPLICIT NONE
     CONTAINS

     ! THIS WORK WAS DESCRIBED IN DETAIL IN: S. A. TRYGUBENKO AND D. J. WALES, `A DOUBLY NUDGED ELASTIC BAND METHOD FOR FINDING
     ! TRANSITION STATES', J. CHEM. PHYS., 120, 2082-2094 (2004). SUMMARY IS AVAILABLE ONLINE AT
     ! HTTP://WWW-WALES.CH.CAM.AC.UK/~SAT39/DNEBTESTS/ 
     SUBROUTINE NEBGRADIENT
          USE KEYGRAD
          USE KEYTAU
          USE NEBDATA
          USE KEYNEB,ONLY: NIMAGE
          USE TANGENT
          USE SPFUNCTS
          USE NEBUTILS
          USE MODUNRES, ONLY: C,NRES
          USE KEY, ONLY: UNRST, FROZEN, FREEZE, NEBK, BULKT, NEBRESEEDT, NEBRESEEDEMAX, NEBRESEEDBMAX, &
  &                      RBAAT, NREPULSIVE, NEBRESEEDINT, &
  &                      ORDERI,ORDERJ,EPSALPHA, REPPOW,PERMDIST, DEBUG, TWOD, RIGIDBODY, NEBRESEEDDEL1, &
  &                      DISTREF, ADDREPT, NEBRESEEDDEL2, NEBRESEEDPOW1, NEBRESEEDPOW2
          USE COMMONS, ONLY : PARAM1, PARAM2, PARAM3
          IMPLICIT NONE
           
          INTEGER :: J1,J2,K,J3,J4,J5,JOUT,JIN,NSHRINK,NSHRINKATOM(NATOMS),J6,NATTACH,NATTACHS
          DOUBLE PRECISION :: PERPCOMP=0.1D0
          DOUBLE PRECISION :: DIHEDIST,TMPRMS,QINT(NINTS*(NIMAGE+2)),DIFFM(NINTS),DIFFP(NINTS) ! JMC
          DOUBLE PRECISION DIST
          DOUBLE PRECISION :: PI=3.141592653589793D0
          DOUBLE PRECISION TEMP1(NOPT), TEMP2(NOPT), GGGSAVE(NOPT*(NIMAGE+2)), DUMMY, DUMMY2, SPRING(NOPT*(NIMAGE+2)), &
  &                        TRUEPERP(NOPT*(NIMAGE+2)), IMCOORDS(3*NATOMS), RSITE(3*NATOMS), RMAT(3,3), D, DIST2, LX(3), LV(3), &
  &                        REPGRAD(3), MAGREP, MEAND, VS(3), VF(3), MP1S(3), MP2S(3), MP1F(3), MP2F(3), DS, DF, D1, D2, DI, DJ, &
  &                        MP2SS(3), MP2FF(3), MP1SS(3), MP1FF(3), DSS, DFF, DISTA, DISTB
          INTEGER, ALLOCATABLE :: IREPTEMP(:)
          LOGICAL SHRINKT(NATOMS)
          
          CALL TRUEPOTEG(.TRUE.)

          CALL MEPTANGENT

          IF (BADTAU) RETURN

!
! RBAAT ACTUALLY WORKS BEST WITH STANDARD DNEB PROCEDURE. HMM.
!

!
! FOR RBAAT RIGID BODIES FREEZE ONE SET OF COORDINATES BY SETTING GRADIENT
! COMPONENTS TO ZERO. WE ARE THEN JUST MINIMISING WITH RESPECT TO THE OTHER SET
! FOR EACH IMAGE.
!
!         IF (RBAAT) THEN
!            DO J1=2,NIMAGE+1
!               GGG(NOPT*(J1-1)+1:NOPT*(J1-1)+(NOPT/2))=0.0D0 ! FREEZE CENTRES OF MASS
!               GGG(NOPT*(J1-1)+(NOPT/2)+1:NOPT*J1)=0.0D0 ! FREEZE ANGLES
!            ENDDO
!            GOTO 555
!         ENDIF
!          
! GGGSAVE CONTAINS THE TRUE GRADIENT ON IMAGES
!
          GGGSAVE(1:NOPT*(NIMAGE+2))=GGG(1:NOPT*(NIMAGE+2))
          ! GGG CONTAINS THE PERPENDICULAR COMPONENT OF THE TRUE GRADIENT AFTER THIS BLOCK
          DO J1=2,NIMAGE+1  !  MEP-PERPENDICULAR COMPONENT OF TRUE GRADIENT ONLY: G = G - G|| = G - (G,TAU)*TAU
               IF (UNRST) THEN
                  GGG(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS) = GGG(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS) - &
               &    DOT_PRODUCT(GGG(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS),TANVEC(:NINTS,J1-1))*TANVEC(:NINTS,J1-1)
               ELSE
                  GGG(NOPT*(J1-1)+1:NOPT*J1) = GGG(NOPT*(J1-1)+1:NOPT*J1) - &
               &    DOT_PRODUCT(GGG(NOPT*(J1-1)+1:NOPT*J1),TANVEC(:,J1-1))*TANVEC(:,J1-1)
               ENDIF
               TRUEPERP(1:NOPT*(NIMAGE+2))=GGG(1:NOPT*(NIMAGE+2))
          ENDDO
          IF (UNRST) THEN
             TMPRMS=0.0D0
             DO J1=1,NIMAGE
                TMPRMS = TMPRMS + DOT_PRODUCT(G(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS),G(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS))
             ENDDO
             RMS=SQRT( TMPRMS/(NIMAGE*NINTS) ) 
          ELSE
             RMS=SQRT( DOT_PRODUCT(G,G)/(NIMAGE*NOPT) ) 
          ENDIF

          IF (BULKT.AND.(GRADTYPE.NE.'DNEB')) THEN
             PRINT '(A)',' NEBGRADIENT> ERROR - MINIMUM IMAGE DISTANCES ONLY CODED FOR GRADTYPE DNEB'
             STOP
          ENDIF
          
          ! SPRING GRADIENT ON IMAGES
          SELECT CASE (GRADTYPE) ! THIS CALCULATES SPRING GRADIENT
             CASE("JNEW") ! SPRINGS IN NEW IMPLEMENTATION OF NEB BY JONSSON
               DO J1=1,NIMAGE
                 GSPR(NOPT*(J1-1)+1:NOPT*J1) = &
                 & - NEBK*( SQRT(SUM( ( XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) - XYZ(NOPT*J1+1:NOPT*(J1+1)) )**2 )) &
                 &          - SQRT(SUM( ( XYZ(NOPT*(J1-1)+1:NOPT*J1)     - XYZ(NOPT*J1+1:NOPT*(J1+1)) )**2 )) )*TANVEC(:,J1)
               ENDDO
             CASE("SPR") ! FULL SPRINGS:  GSPR = K [2X(J) -X(J-1) -X(J+1)]
               DO J1=1,NIMAGE
                 GSPR(NOPT*(J1-1)+1:NOPT*J1) = &
                 & NEBK*( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) )
               ENDDO
             CASE("JOLD") ! SPRINGS MEP-PARALLEL COMPONENT ONLY (ORIGINAL FORMULATION): G = (GSPR,TAU)*TAU
               DO J1=1,NIMAGE
                 GSPR(NOPT*(J1-1)+1:NOPT*J1) = &
                 & NEBK*DOT_PRODUCT( (2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2))), &
                 & TANVEC(:,J1) )*TANVEC(:,J1)
               ENDDO
              CASE("NPER")! SPRINGS MEP-PARALLEL COMPONENT (ORIGINAL FORMULATION) + N% OF PERPENDICULAR COMPONENT:
               DO J1=1,NIMAGE   !GJ = GJ + G|_*0.1 + GSPR|| == GJ + [ GSPR - (GSPR,TAU)*TAU ]*N + ...
                 GSPR(NOPT*(J1-1)+1:NOPT*J1) = ( &
                 &   NEBK*           ( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) )& 
                 & - NEBK*DOT_PRODUCT( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) ,&
                 & TANVEC(:,J1) )*TANVEC(:,J1) &
                                               & )*PERPCOMP
                 
                 GSPR(NOPT*(J1-1)+1:NOPT*J1) = GSPR(NOPT*(J1-1)+1:NOPT*J1) &  !... (GSPR,TAU)*TAU
                 & + NEBK*DOT_PRODUCT( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) ,&
                 & TANVEC(:,J1) )*TANVEC(:,J1)
               ENDDO
              CASE("NPER2")! SAME AS ABOVE BUT OPTIMISED FOR EFFICIENCY?
               DO J1=1,NIMAGE
                 ! CALCULATE GS (SPRING GRADIENT ON IMAGE J1)
                 GSPR(NOPT*(J1-1)+1:NOPT*J1) = NEBK*(&
                 2*XYZ(NOPT*J1+1:NOPT*(J1+1))-XYZ(NOPT*(J1-1)+1:NOPT*J1)-XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) &
                                                       &)
                 ! CALCULATE [ F*GS + (1-F)(GS,TAU)TAU ]
                 GSPR(NOPT*(J1-1)+1:NOPT*J1) = &
                 & PERPCOMP*GSPR(NOPT*(J1-1)+1:NOPT*J1) &
                 & +(1-PERPCOMP)*DOT_PRODUCT(GSPR(NOPT*(J1-1)+1:NOPT*J1),TANVEC(:,J1))*TANVEC(:,J1)
               ENDDO
              CASE("DNEB")  ! THIS APPEARS TO BE THE DEFAULT
               DO J1=1,NIMAGE
                 ! CALCULATE SSS = ~G PERP IN DNEB PAPER, THE PERPENDICULAR PART OF THE SPRING GRADIENT VECTOR
!
!  HERE WE ARE CALCULATING SPRING GRADIENT VECTOR - (SPRING GRADIENT VECTOR . TANGENT VECTOR) TANGENT VECTOR
!  THIS GIVES THE SPRING GRADIENT PERPENDICULAR COMPONENT, WHICH IS PUT IN SSS.
!
!  CHANGED TO USE DYNAMICALLY ADJUSTED NEBK DJW 21/10/08
!  ONLY DONE FOR DEFAULT CASE("DNEB")
!
                  IF (BULKT) THEN ! MINIMUM IMAGE CONVENTION FOR DISTANCES!
                     DO K=1,NATOMS
                        TEMP1(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1 ) &
   &                       -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1))/PARAM1)
                        TEMP1(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2 ) &
   &                       -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2))/PARAM2)
                        IF (.NOT.TWOD) TEMP1(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3 ) &
   &                       -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3))/PARAM3)

                        TEMP2(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1+1)+3*(K-1)+1 ) &
   &                       -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1+1)+3*(K-1)+1))/PARAM1)
                        TEMP2(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1+1)+3*(K-1)+2 ) &
   &                       -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1+1)+3*(K-1)+2))/PARAM2)
                        IF (.NOT.TWOD) TEMP2(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1+1)+3*(K-1)+3 ) &
   &                       -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1+1)+3*(K-1)+3))/PARAM3)
                     ENDDO
                     SSS(NOPT*J1+1:NOPT*(J1+1)) = NEWNEBK(J1)*TEMP1(1:NOPT)+NEWNEBK(J1+1)*TEMP2(1:NOPT) &
   &                    -NEWNEBK(J1)  *DOT_PRODUCT(TEMP1(1:NOPT),TANVEC(1:NOPT,J1) )*TANVEC(1:NOPT,J1) &
   &                    -NEWNEBK(J1+1)*DOT_PRODUCT(TEMP2(1:NOPT),TANVEC(1:NOPT,J1) )*TANVEC(1:NOPT,J1)
                  ELSE
                     SPRING(NOPT*J1+1:NOPT*(J1+1))=  &
      &   NEWNEBK(J1)  *( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) ) & 
      & + NEWNEBK(J1+1)*( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) ) 

                     SSS(NOPT*J1+1:NOPT*(J1+1)) = &
      &   NEWNEBK(J1)  *( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) ) &
      & + NEWNEBK(J1+1)*( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) ) &
      & - NEWNEBK(J1)  *DOT_PRODUCT( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1)    , TANVEC(:,J1) )*TANVEC(:,J1) &
      & - NEWNEBK(J1+1)*DOT_PRODUCT( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)), TANVEC(:,J1) )*TANVEC(:,J1)
                  ENDIF
!                 PRINT '(A,I6)','GRAD> SSS SPRING GRADIENT PERPENDICULAR FOR IMAGE ',J1
!                 PRINT '(6G20.10)',SSS(1+NOPT*J1:(J1+1)*NOPT)

!
! SSS NOW CONTAINS THE SPRING GRADIENT PERPENDICULAR PART.
! SPRING CONTAINS THE COMPLETE SPRING GRADIENT.
!
                  CALL NUDGE(J1,J1)
!                 PRINT '(A,I6)','GRAD> SSS  G~* FOR IMAGE ',J1
!                 PRINT '(6G20.10)',SSS(1+NOPT*J1:(J1+1)*NOPT)
!
! ON EXIT FROM NUDGE SSS CONTAINS G~* FROM EQUATION (13) OF THE DNEB PAPER.
! GGG STILL CONTAINS G PERP, THE PERPENDICULAR PART OF THE TRUE GRADIENT.
!
               ENDDO
!
! THE CALL TO NUDGE CALCULATED A NEW SSS CONTAINING THE NUDGED SPRING GRADIENT G~*.
! WE NOW ADD THIS TO GGG. WE ARE STILL MISSING ~G PARALLEL.
!
               GGG(NOPT+1:NOPT*(NIMAGE+1)) = GGG(NOPT+1:NOPT*(NIMAGE+1)) + SSS(NOPT+1:NOPT*(NIMAGE+1))
!                 PRINT '(A,I6)','GRAD> GGG G PERP + G~* FOR IMAGE ',J1
!                 PRINT '(6G20.10)',GGG(1+NOPT*J1:(J1+1)*NOPT)
!
! CALCULATE ~G PARALLEL AND STORE IN SSS. THIS IS NOT ADDED TO GGG UNTIL AFTER THE WHOLE SELECT BLOCK.
! IN FACT ~G PARALLEL IS NOT USED AS THE FORMULA EXPECTED FROM PROJECTION, BUT USES AN
! ALTERNATIVE FORMULATION FROM EQUATION (5) OF THE DNEB PAPER.
!
               DO J1=1,NIMAGE 
                  IF (BULKT) THEN ! MINIMUM IMAGE CONVENTION FOR DISTANCES!

                     DISTA=0.0D0; DISTB=0.0D0
                     DO J2=1,NATOMS
                        DISTA=DISTA+ &
  &                           (XYZ(NOPT*(J1+1)+3*(J2-1)+1)-XYZ(NOPT*J1+3*(J2-1)+1) - &
  &               PARAM1*NINT((XYZ(NOPT*(J1+1)+3*(J2-1)+1)-XYZ(NOPT*J1+3*(J2-1)+1))/PARAM1))**2+ &
  &                           (XYZ(NOPT*(J1+1)+3*(J2-1)+2)-XYZ(NOPT*J1+3*(J2-1)+2) - &
  &               PARAM2*NINT((XYZ(NOPT*(J1+1)+3*(J2-1)+2)-XYZ(NOPT*J1+3*(J2-1)+2))/PARAM2))**2
   IF (.NOT.TWOD) DISTA=DISTA+(XYZ(NOPT*(J1+1)+3*(J2-1)+3)-XYZ(NOPT*J1+3*(J2-1)+3) - &
  &               PARAM3*NINT((XYZ(NOPT*(J1+1)+3*(J2-1)+3)-XYZ(NOPT*J1+3*(J2-1)+3))/PARAM3))**2
                        DISTB=DISTB+ &
  &                           (XYZ(NOPT*(J1-1)+3*(J2-1)+1)-XYZ(NOPT*J1+3*(J2-1)+1) - &
  &               PARAM1*NINT((XYZ(NOPT*(J1-1)+3*(J2-1)+1)-XYZ(NOPT*J1+3*(J2-1)+1))/PARAM1))**2+ &
  &                           (XYZ(NOPT*(J1-1)+3*(J2-1)+2)-XYZ(NOPT*J1+3*(J2-1)+2) - &
  &               PARAM2*NINT((XYZ(NOPT*(J1-1)+3*(J2-1)+2)-XYZ(NOPT*J1+3*(J2-1)+2))/PARAM2))**2
   IF (.NOT.TWOD) DISTA=DISTA+(XYZ(NOPT*(J1-1)+3*(J2-1)+3)-XYZ(NOPT*J1+3*(J2-1)+3) - &
  &               PARAM3*NINT((XYZ(NOPT*(J1-1)+3*(J2-1)+3)-XYZ(NOPT*J1+3*(J2-1)+3))/PARAM3))**2
                     ENDDO
                     SSS(NOPT*J1+1:NOPT*(J1+1))= - ( NEWNEBK(J1+1)*SQRT(DISTA) - NEWNEBK(J1)*SQRT(DISTB) )*TANVEC(1:NOPT,J1)
                  ELSE
                     SSS(NOPT*J1+1:NOPT*(J1+1))= &
  &               - ( NEWNEBK(J1+1)*SQRT(SUM( ( XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) - XYZ(NOPT*J1+1:NOPT*(J1+1)) )**2 )) &
  &                 - NEWNEBK(J1)  *SQRT(SUM( ( XYZ(NOPT*(J1-1)+1:NOPT*J1)     - XYZ(NOPT*J1+1:NOPT*(J1+1)) )**2 )) )*TANVEC(:,J1)
                  ENDIF
               ENDDO
              CASE("DNEB2") ! SAME AS "DNEB", EXCEPT SPRING GRADIENT USAGE IS MORE CONSISTENT;
               DO J1=1,NIMAGE ! THANKS TO DR.~DOMINIC R. ALFONSO FOR POINTING THAT OUT
                 ! CALCULATE GSPR PERP
                 SSS(NOPT*J1+1:NOPT*(J1+1)) = &
                 &   NEBK*           ( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) )&
                 & - NEBK*DOT_PRODUCT( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) ,&
                 & TANVEC(:,J1) )*TANVEC(:,J1)
                 CALL NUDGE(J1,J1)
               ENDDO
!
! GGG CONTAINS THE TRUE GRADIENT COMPONENT PERPENDICULAR AND SSS NOW CONTAINS G~*,
! SO THE NEXT LINE GIVES THE PERPENDICULAR COMPONENT OF THE DNEB GRADIENT.
!
               GGG(NOPT+1:NOPT*(NIMAGE+1)) = GGG(NOPT+1:NOPT*(NIMAGE+1)) + SSS(NOPT+1:NOPT*(NIMAGE+1))

               DO J1=1,NIMAGE ! NOW CALCULATE GSPR PARALLEL (USING OLD FORMULA)
                  IF (BULKT) THEN ! MINIMUM IMAGE CONVENTION FOR DISTANCES!
                     DO K=1,NATOMS
                        TEMP1(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1 ) &
   &                       -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1))/PARAM1)
                        TEMP1(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2 ) &
   &                       -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2))/PARAM2)
                        IF (.NOT.TWOD) TEMP1(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3 ) &
   &                       -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3))/PARAM3)

                        TEMP2(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1+1)+3*(K-1)+1 ) &
   &                       -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1+1)+3*(K-1)+1))/PARAM1)
                        TEMP2(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1+1)+3*(K-1)+2 ) &
   &                       -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1+1)+3*(K-1)+2))/PARAM2)
                        IF (.NOT.TWOD) TEMP2(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1+1)+3*(K-1)+3 ) &
   &                       -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1+1)+3*(K-1)+3))/PARAM3)
                     ENDDO
                     SSS(NOPT*J1+1:NOPT*(J1+1)) = &
   &                    NEBK*DOT_PRODUCT(TEMP1(1:NOPT)+TEMP2(1:NOPT),TANVEC(1:NOPT,J1) )*TANVEC(1:NOPT,J1)
                  ELSE
                     SSS(NOPT*J1+1:NOPT*(J1+1)) = &
                 & NEBK*DOT_PRODUCT( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) ,&
                 &   TANVEC(:,J1) )*TANVEC(:,J1)
                  ENDIF
               ENDDO
              CASE("DIHEN") ! DAE DISTANCE BETWEEN IMAGES CALCULATED IN DIHEDRAL ANGLE SPACE; FOR USE WITH CHARMM AND UNRES.
                            ! USE NEW FORMULATION OF SPRING FORCE TO SIMPLIFY
               SSS=0.0D0    ! INITIALISE SPRING GRADIENT
               DO J1=1,NIMAGE
                 IF (UNRST) THEN
                    CALL UNRESGETDIHEDIST(DIHEDIST,XYZ(NOPT*(J1-1)+1:NOPT*J1),XYZ(NOPT*J1+1:NOPT*(J1+1)),&
                    &XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)))
                 ELSE
                    CALL GETDIHEDIST(DIHEDIST,XYZ(NOPT*(J1-1)+1:NOPT*J1),XYZ(NOPT*J1+1:NOPT*(J1+1)),&
                    &XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)))
                 ENDIF
                 GSPR(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS) = - NEBK * DIHEDIST * TANVEC(:NINTS,J1)
               ENDDO
              CASE("DNEBU") ! DNEB FOR UNRES
               SSS=0.0D0    ! INITIALISE SPRING GRADIENT
               DO J1=1,NIMAGE+2 ! GET INTERNALS FOR ENDPOINTS AND IMAGES
                  DO J2=1,NRES
                     C(1,J2)=XYZ(6*(J2-1)+1+NOPT*(J1-1))
                     C(2,J2)=XYZ(6*(J2-1)+2+NOPT*(J1-1))
                     C(3,J2)=XYZ(6*(J2-1)+3+NOPT*(J1-1))
                     C(1,J2+NRES)=XYZ(6*(J2-1)+4+NOPT*(J1-1))
                     C(2,J2+NRES)=XYZ(6*(J2-1)+5+NOPT*(J1-1))
                     C(3,J2+NRES)=XYZ(6*(J2-1)+6+NOPT*(J1-1))
                  ENDDO
                  CALL UPDATEDC
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL GEOM_TO_VAR(NINTS,QINT(NINTS*(J1-1)+1:NINTS*J1))
               ENDDO

               DO J1=1,NIMAGE
                 ! CALCULATE GSPR PERP
                 DIFFM(1:NINTS) = QINT(NINTS*J1+1:NINTS*(J1+1)) - QINT(NINTS*(J1-1)+1:NINTS*J1)
                 DIFFP(1:NINTS) = QINT(NINTS*J1+1:NINTS*(J1+1)) - QINT(NINTS*(J1+1)+1:NINTS*(J1+2))
                 DO J2=1,NINTS
                    IF (DIFFM(J2).LT.-PI) DIFFM(J2) = DIFFM(J2)+2.0D0*PI
                    IF (DIFFP(J2).GT.PI) DIFFP(J2) = DIFFP(J2)-2.0D0*PI
                 ENDDO

                 SSS(NOPT*J1+1:NOPT*J1+NINTS) = NEBK*( DIFFM(1:NINTS) + DIFFP(1:NINTS) ) &
                 & - NEBK*DOT_PRODUCT( ( DIFFM(1:NINTS) + DIFFP(1:NINTS) ) , TANVEC(1:NINTS,J1) )*TANVEC(1:NINTS,J1)
                 CALL NUDGE(J1,J1) ! ZERO-ORDER NUDGING
               ENDDO

               GGG(NOPT+1:NOPT*(NIMAGE+1)) = GGG(NOPT+1:NOPT*(NIMAGE+1)) + SSS(NOPT+1:NOPT*(NIMAGE+1)) ! G PERP + G TILDE STAR

               DO J1=1,NIMAGE ! NOW CALCULATE GSPR PARALLEL AS FOR DIHEN
                 IF (UNRST) THEN
                    CALL UNRESGETDIHEDIST(DIHEDIST,XYZ(NOPT*(J1-1)+1:NOPT*J1),XYZ(NOPT*J1+1:NOPT*(J1+1)),&
                    &XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)))
                 ELSE
                    CALL GETDIHEDIST(DIHEDIST,XYZ(NOPT*(J1-1)+1:NOPT*J1),XYZ(NOPT*J1+1:NOPT*(J1+1)),&
                    &XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)))
                 ENDIF
                 SSS(NOPT*J1+1:NOPT*J1+NINTS) = - NEBK * DIHEDIST * TANVEC(:NINTS,J1)
                 !PRINT *,'GSPR ',SSS(NOPT*J1+1:NOPT*J1+NINTS)
               ENDDO
              CASE("DNEB3")  ! USES DEVIATION FROM THE AVERAGE SEPARATION IN THE SPRING TERM INSTEAD OF ABSOLUTE DISTANCE
               CALL DISTANCES ! SETS UP IMAGE DISTANCES IN VECTOR DVEC FROM 1 TO NIMAGE+1
               MEAND=0.0D0
               DO J1=1,NIMAGE+1
                  MEAND=MEAND+DVEC(J1)
               ENDDO
               MEAND=MEAND/(NIMAGE+1)
               PRINT '(A,F15.3)','GRAD> MEAN IMAGE SEPARATION=',MEAND
               
               DO J1=1,NIMAGE
                 ! CALCULATE SSS = GSPR PERP  I.E. ~G PERP IN DNEB PAPER
!
!  HERE WE ARE CALCULATING SPRING GRADIENT VECTOR - (SPRING GRADIENT VECTOR . TANGENT VECTOR) TANGENT VECTOR
!  THIS GIVES THE SPRING GRADIENT PERPENDICULAR COMPONENT
!
                  IF (BULKT) THEN ! MINIMUM IMAGE CONVENTION FOR DISTANCES!
                     DO K=1,NATOMS
                        TEMP1(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1 ) &
   &                       -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1))/PARAM1)
                        TEMP1(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2 ) &
   &                       -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2))/PARAM2)
                        IF (.NOT.TWOD) TEMP1(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3 ) &
   &                       -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3))/PARAM3)

                        TEMP2(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1+1)+3*(K-1)+1 ) &
   &                       -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1+1)+3*(K-1)+1))/PARAM1)
                        TEMP2(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1+1)+3*(K-1)+2 ) &
   &                       -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1+1)+3*(K-1)+2))/PARAM2)
                        IF (.NOT.TWOD) TEMP2(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1+1)+3*(K-1)+3 ) &
   &                       -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1+1)+3*(K-1)+3))/PARAM3)
                     ENDDO
                     SSS(NOPT*J1+1:NOPT*(J1+1)) = NEWNEBK(J1)*TEMP1(1:NOPT)*(DVEC(J1)-MEAND)/DVEC(J1) &
   &                                             +NEWNEBK(J1+1)*TEMP2(1:NOPT)*(DVEC(J1+1)-MEAND)/DVEC(J1+1) &
   &                   -NEWNEBK(J1)  *DOT_PRODUCT(TEMP1(1:NOPT)*(DVEC(J1)-MEAND)/DVEC(J1),TANVEC(1:NOPT,J1) )*TANVEC(1:NOPT,J1) &
   &                   -NEWNEBK(J1+1)*DOT_PRODUCT(TEMP2(1:NOPT)*(DVEC(J1+1)-MEAND)/DVEC(J1+1),TANVEC(1:NOPT,J1) )*TANVEC(1:NOPT,J1)
                  ELSE
                     SSS(NOPT*J1+1:NOPT*(J1+1)) = &
      &   NEWNEBK(J1)  *( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) )*(DVEC(J1)-MEAND)/DVEC(J1) &
      & + NEWNEBK(J1+1)*( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) )*(DVEC(J1+1)-MEAND)/DVEC(J1+1) &
      & - NEWNEBK(J1)  *DOT_PRODUCT( (XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1))*(DVEC(J1)-MEAND)/DVEC(J1)    &
                                    , TANVEC(:,J1) )*TANVEC(:,J1) &
      & - NEWNEBK(J1+1)*DOT_PRODUCT( (XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)))*(DVEC(J1+1)-MEAND)/DVEC(J1+1) &
      &                             , TANVEC(:,J1) )*TANVEC(:,J1)
                  ENDIF
                  SPRING(NOPT*(J1-1)+1:NOPT*J1)=SSS(NOPT*(J1-1)+1:NOPT*J1)
                  CALL NUDGE(J1,J1)
               ENDDO
!
! HERE WE COMPLETE THE DNEB GRADIENT VECTOR AS IN EQUATION (12) OF THE DNEB PAPER.
! THE CALL TO NUDGE CALCULATED A NEW SSS CONTAINING THE NUDGED SPRING GRADIENT.
!
               GGG(NOPT+1:NOPT*(NIMAGE+1)) = GGG(NOPT+1:NOPT*(NIMAGE+1)) + SSS(NOPT+1:NOPT*(NIMAGE+1))

               DO J1=1,NIMAGE ! NOW CALCULATE GSPR PARALLEL (USING FORMULA THAT KEEPS IMAGES EQUISPACED WITH NEW TANGENT)
                 SSS(NOPT*J1+1:NOPT*(J1+1)) = &
                 & - ( NEWNEBK(J1+1)*SQRT(SUM( ( XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) - XYZ(NOPT*J1+1:NOPT*(J1+1)) )**2 )) &
                 &   - NEWNEBK(J1)  *SQRT(SUM( ( XYZ(NOPT*(J1-1)+1:NOPT*J1)     - XYZ(NOPT*J1+1:NOPT*(J1+1)) )**2 )) )*TANVEC(:,J1)
               ENDDO
              CASE DEFAULT
               PRINT *,'ERROR: UNKNOWN GRADIENT TYPE "'//TRIM(ADJUSTL(GRADTYPE))//'"'; STOP
          END SELECT

!----------------------------------------------------------------------------------------------------------------------------------
!
! END OF LARGE SELECT BLOCK. NOW WE ADD ON SSS TO GGG, WHICH COMPLETES EQUATION (12) OF THE DNEB PAPER.
! AT THIS POINT GG CONTAINS G PERP + ~G* AND WE NEED TO ADD ON ~G PARALLEL.
!

          GGG(NOPT+1:NOPT*(NIMAGE+1)) = GGG(NOPT+1:NOPT*(NIMAGE+1)) + SSS(NOPT+1:NOPT*(NIMAGE+1))

          IF (NEBRESEEDT.AND.(NREPULSIVE.GT.0).AND.ADDREPT) THEN
             DO J2=1,NREPULSIVE
                DO J1=2,NIMAGE+1
                   IF ((BADIMAGE(J1).AND.(REPPOW(J2).GT.0)).OR.(BADPEPTIDE(J1).AND.(REPPOW(J2).LT.0))) THEN
!
!  ADD REPULSIVE/CONSTRAINT TERMS DEFINED AS INTERATOMIC DISTANCES.
!

                      DIST=SQRT((XYZ(NOPT*(J1-1)+3*(ORDERI(J2)-1)+1)-XYZ(NOPT*(J1-1)+3*(ORDERJ(J2)-1)+1))**2 &
  &                            +(XYZ(NOPT*(J1-1)+3*(ORDERI(J2)-1)+2)-XYZ(NOPT*(J1-1)+3*(ORDERJ(J2)-1)+2))**2 &
  &                            +(XYZ(NOPT*(J1-1)+3*(ORDERI(J2)-1)+3)-XYZ(NOPT*(J1-1)+3*(ORDERJ(J2)-1)+3))**2)

                      REPGRAD(1:3*NATOMS)=0.0D0

                      DUMMY=DIST-DISTREF(J2)  !  DJW BUG THIS J2 WAS J1!
                      IF (DUMMY.EQ.0.0D0) DUMMY=1.0D-10

                      REPGRAD(3*(ORDERI(J2)-1)+1:3*(ORDERI(J2)-1)+3)= &
     &                 -REPPOW(J2)*EPSALPHA(J2) &
     &                  *(XYZ((J1-1)*NOPT+3*(ORDERI(J2)-1)+1:(J1-1)*NOPT+3*(ORDERI(J2)-1)+3) &
     &                   -XYZ((J1-1)*NOPT+3*(ORDERJ(J2)-1)+1:(J1-1)*NOPT+3*(ORDERJ(J2)-1)+3)) &
     &                         /(DIST*DUMMY**(REPPOW(J2)+1))

                      REPGRAD(3*(ORDERJ(J2)-1)+1:3*(ORDERJ(J2)-1)+3)= &
     &                 -REPPOW(J2)*EPSALPHA(J2) &
     &                  *(XYZ((J1-1)*NOPT+3*(ORDERJ(J2)-1)+1:(J1-1)*NOPT+3*(ORDERJ(J2)-1)+3) &
     &                   -XYZ((J1-1)*NOPT+3*(ORDERI(J2)-1)+1:(J1-1)*NOPT+3*(ORDERI(J2)-1)+3)) &
     &                         /(DIST*DUMMY**(REPPOW(J2)+1))

                      GGG((J1-1)*NOPT+3*(ORDERI(J2)-1)+1:(J1-1)*NOPT+3*(ORDERI(J2)-1)+3)= &
    &                 GGG((J1-1)*NOPT+3*(ORDERI(J2)-1)+1:(J1-1)*NOPT+3*(ORDERI(J2)-1)+3)  &
    &                   +REPGRAD(3*(ORDERI(J2)-1)+1:3*(ORDERI(J2)-1)+3)
  
                      GGG((J1-1)*NOPT+3*(ORDERJ(J2)-1)+1:(J1-1)*NOPT+3*(ORDERJ(J2)-1)+3)= &
    &                 GGG((J1-1)*NOPT+3*(ORDERJ(J2)-1)+1:(J1-1)*NOPT+3*(ORDERJ(J2)-1)+3)  &
    &                   +REPGRAD(3*(ORDERJ(J2)-1)+1:3*(ORDERJ(J2)-1)+3)

                   ENDIF
                ENDDO
             ENDDO
          ENDIF
!
! SET GRADIENTS ON FROZEN ATOMS TO ZERO.
!
555       IF (FREEZE) THEN
             DO J1=2,NIMAGE+1  
                DO J2=1,NATOMS
                   IF (FROZEN(J2)) THEN
                      GGG(NOPT*(J1-1)+3*(J2-1)+1)=0.0D0
                      GGG(NOPT*(J1-1)+3*(J2-1)+2)=0.0D0
                      GGG(NOPT*(J1-1)+3*(J2-1)+3)=0.0D0
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          RMS=0.0D0
          DO J1=2,NIMAGE+1
             DO J2=1,NOPT
                RMS=RMS+GGG(NOPT*(J1-1)+J2)**2
             ENDDO
          ENDDO
          RMS=SQRT(RMS/(NOPT*NIMAGE))

          ETOTAL=SUM(EIMAGE)

          IF (ORT) THEN ! REMOVE OVERALL ROTATION/TRANSLATION BY FREEZING 6 DEGREES OF FREEDOM
               G(1::NOPT)=0.0D0
               G(2::NOPT)=0.0D0
               G(3::NOPT)=0.0D0
               G(4::NOPT)=0.0D0
               G(5::NOPT)=0.0D0
               G(7::NOPT)=0.0D0
          ENDIF
     END SUBROUTINE NEBGRADIENT
          
     SUBROUTINE NUDGE(JS,JT) ! SUBROUTINE DOES HIGHER-ORDER NUDGINGS; IT MODIFIES GSPR, EATS BOTH GSPR AND G; PREREQ.
          ! - GSPR AND G ARE CALCULATED; JS - IMAGE ON WHICH GSPER IS NUDGED BY GTPER ON IMAGE JT
          USE NEBDATA
          USE KEY, ONLY: UNRST
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: JS,JT
!
! ON ENTRY, SSS CONTAINS THE PERPENDICULAR SPRING GRADIENT COMPONENT
! GGG MUST CONTAIN A NON-UNIT VECTOR WITH THE PERPENDICULAR COMPONENT OF THE TRUE POTENTIAL
! ON EXIT, SSS CONTAINS G~* FROM EQUATION (13) OF THE DNEB PAPER
!
          IF (UNRST) THEN
             IF (DOT_PRODUCT(SSS(NOPT*JS+1:NOPT*JS+NINTS),GGG(NOPT*JT+1:NOPT*JT+NINTS)) < 0.0D0) THEN
                  ! CALCULATE THE BEAST GSPR PERP - (GSPR PERP, GT PERP)* GT PERP/(GT PERP, GT PERP)
                  SSS(NOPT*JS+1:NOPT*JS+NINTS) = SSS(NOPT*JS+1:NOPT*JS+NINTS) - &
                  & DOT_PRODUCT(SSS(NOPT*JS+1:NOPT*JS+NINTS),GGG(NOPT*JT+1:NOPT*JT+NINTS))*GGG(NOPT*JT+1:NOPT*JT+NINTS)/ &
                  & DOT_PRODUCT(GGG(NOPT*JT+1:NOPT*JT+NINTS),GGG(NOPT*JT+1:NOPT*JT+NINTS))
             ENDIF
          ELSE
             IF (DOT_PRODUCT(SSS(NOPT*JS+1:NOPT*(JS+1)),GGG(NOPT*JT+1:NOPT*(JT+1))) < 0.0D0) THEN
                  ! CALCULATE THE BEAST GSPR PERP - (GSPR PERP, GT PERP)* GT PERP/(GT PERP, GT PERP)
                  SSS(NOPT*JS+1:NOPT*(JS+1)) = SSS(NOPT*JS+1:NOPT*(JS+1)) - &
                  & DOT_PRODUCT(SSS(NOPT*JS+1:NOPT*(JS+1)),GGG(NOPT*JT+1:NOPT*(JT+1)))*GGG(NOPT*JT+1:NOPT*(JT+1))/ &
                  & DOT_PRODUCT(GGG(NOPT*JT+1:NOPT*(JT+1)),GGG(NOPT*JT+1:NOPT*(JT+1)))
             ENDIF
          ENDIF
     END SUBROUTINE NUDGE

     SUBROUTINE TRUEPOTEG(GTEST)
       USE PORFUNCS
       USE NEBDATA
       USE KEYNEB,ONLY: NIMAGE
       USE MODUNRES, ONLY: C,NRES
       USE KEY, ONLY: UNRST, FREEZENODEST, INTEPSILON
       USE INTCOMMONS, ONLY : DESMINT, NNZ, KD, NINTC, PREVDIH, DIHINFO
       USE MODCHARMM, ONLY: CHRMMT, ACESOLV, NCHENCALLS, ACEUPSTEP
       IMPLICIT NONE

       LOGICAL, INTENT(IN) :: GTEST

       INTEGER :: I,J
       DOUBLE PRECISION :: RMSTMP,HARVEST,DPRAND

       DOUBLE PRECISION :: DUMINT(NINTC)

       ! ENERGY AND GRADIENT FOR IMAGES
!      PRINT '(A)',' '
!      PRINT '(A,I8,A,G20.10,A)','IMAGE ',1,' ENERGY ',EEE(1),' POINTS:'
!      PRINT '(3G20.10)',XYZ(3*NATOMS*(1-1)+1:3*NATOMS*(1-1)+3)
       DO I=2,NIMAGE+1

          IF (FREEZENODEST) THEN
             IF (DESMINT) THEN
                PRINT*, 'NEBGRAD>> DESMINT NOT IMPLEMENTED WITH FREEZENODES.'
                STOP
             ENDIF
             ! FOR NODES THAT DIDN'T MOVE, SET TRUE GRADIENT TO PREVIOUS VALUE
             IF (IMGFREEZE(I-1)) THEN
                GGG(NOPT*(I-1)+1:NOPT*I) = TRUEGRAD(NOPT*(I-1)+1:NOPT*I)
                CYCLE
             ENDIF
          ENDIF

          IF (UNRST) THEN
             IF (DESMINT) THEN
                PRINT*, 'NEBGRAD>> UNRST NOT IMPLEMENTED WITH DESMINT!'
                STOP
             ENDIF
             DO J=1,NRES
                C(1,J)=XYZ(6*(J-1)+1+NOPT*(I-1))
                C(2,J)=XYZ(6*(J-1)+2+NOPT*(I-1))
                C(3,J)=XYZ(6*(J-1)+3+NOPT*(I-1))
                C(1,J+NRES)=XYZ(6*(J-1)+4+NOPT*(I-1))
                C(2,J+NRES)=XYZ(6*(J-1)+5+NOPT*(I-1))
                C(3,J+NRES)=XYZ(6*(J-1)+6+NOPT*(I-1))
             ENDDO
             CALL UPDATEDC
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL CHAINBUILD
          ENDIF

!BS360: UPDATE ACE BORN RADII FOR EACH IMAGE
          IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
          IF (DESMINT) THEN

             CALL POTENTIAL(XYZCART(3*NATOMS*(I-1)+1:3*NATOMS*I),EEE(I), &
  &                   GGGCART(3*NATOMS*(I-1)+1:3*NATOMS*I),.TRUE.,.FALSE.,RMSTMP,.FALSE.,.FALSE.)
             PREVDIH => DIHINFO(I,:)
             CALL TRANSFORM(XYZCART(3*NATOMS*(I-1)+1:3*NATOMS*I),GGGCART(3*NATOMS*(I-1)+1:3*NATOMS*I),&
                  & DUMINT,GGG(NOPT*(I-1)+1:NOPT*I), NINTC,3*NATOMS,NNZ,.TRUE.,.FALSE.,KD,INTEPSILON)
             TRUEGRAD(3*NATOMS*(I-1)+1:3*NATOMS*I)=GGGCART(3*NATOMS*(I-1)+1:3*NATOMS*I) ! SAVE FOR PASSING TO BFGSTS LATER
          ELSE
             CALL POTENTIAL(XYZ(NOPT*(I-1)+1:NOPT*I),EEE(I),GGG(NOPT*(I-1)+1:NOPT*I),.TRUE.,.FALSE.,RMSTMP,.FALSE.,.FALSE.)
             TRUEGRAD(NOPT*(I-1)+1:NOPT*I)=GGG(NOPT*(I-1)+1:NOPT*I) ! SAVE FOR PASSING TO BFGSTS LATER
          ENDIF
!         PRINT '(A,I8,A,G20.10,A)','IMAGE ',I,' ENERGY ',EEE(I),' POINTS:'
!         PRINT '(3G20.10)',XYZ(3*NATOMS*(I-1)+1:3*NATOMS*(I-1)+3)

          IF ( EEE(I) > HUGE(EEE(I)) ) THEN ! BAD GUESS - HIGH-ENERGY IMAGE
             PRINT *, "IMAGE",I," IS BAD! - TRYING TO LOWER IT'S ENERGY..."
             IF (DESMINT) THEN
                PRINT*, 'NEBGRAD>> BAD ENERGY PERTURBATION NOT IMPLEMENTED WITH DESMINT YET.'
                STOP
             ENDIF
             DO J=1,NOPT ! CHANGING GEOMETRY RANDOMLY
                !                        CALL RANDOM_NUMBER(HARVEST)
                HARVEST=DPRAND()
                XYZ(NOPT*(I-1)+J) = XYZ(NOPT*(I-1)+J) + HARVEST*0.01
             ENDDO


             IF (UNRST) THEN
                DO J=1,NRES
                   C(1,J)=XYZ(6*(J-1)+1+NOPT*(I-1))
                   C(2,J)=XYZ(6*(J-1)+2+NOPT*(I-1))
                   C(3,J)=XYZ(6*(J-1)+3+NOPT*(I-1))
                   C(1,J+NRES)=XYZ(6*(J-1)+4+NOPT*(I-1))
                   C(2,J+NRES)=XYZ(6*(J-1)+5+NOPT*(I-1))
                   C(3,J+NRES)=XYZ(6*(J-1)+6+NOPT*(I-1))
                ENDDO
                CALL UPDATEDC
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL CHAINBUILD
             ENDIF

             CALL POTENTIAL(XYZ(NOPT*(I-1)+1:NOPT*I),EEE(I),GGG(NOPT*(I-1)+1:NOPT*I),.TRUE.,.FALSE. &
                  & ,RMSTMP,.FALSE.,.FALSE.)
             TRUEGRAD(NOPT*(I-1)+1:NOPT*I)=GGG(NOPT*(I-1)+1:NOPT*I) ! SAVE FOR PASSING TO BFGSTS LATER

             IF ( EEE(I) > HUGE(EEE(I)) ) THEN
                PRINT *, "FAILED."
                CALL TSUMMARY
                STOP
             ENDIF
          ENDIF
       ENDDO
!      PRINT '(A,I8,A,G20.10,A)','IMAGE ',NIMAGE+2,' ENERGY ',EEE(NIMAGE+2),' POINTS:'
!      PRINT '(3G20.10)',XYZ(3*NATOMS*(NIMAGE+2-1)+1:3*NATOMS*(NIMAGE+2-1)+3)
     END SUBROUTINE TRUEPOTEG

END MODULE GRADIENTS
