      SUBROUTINE TAKESTEP(NP)
      USE commons
      IMPLICIT NONE

      DOUBLE PRECISION DPRAND, RANDOM, XMASS, YMASS, ZMASS, LOCALSTEP, DUMMY2, CDIST(NATOMS), RDOTN, XL, YL, ZL,
     1                 DIST(3*NATOMS), DMAX, VMAX, VMIN, VMAX2, EXPDECAY(NATOMS), XC, YC, ZC, ANGLE, COST, SINT,
     2                 THETA, PHI, PI, DUMMY, CMDIST(NATOMS), CMMAX, RANDOMX, RANDOMY, RANDOMZ, RVEC(3), TX, TY, TZ,
     &                 DELX, DELY, DELZ, SLENGTH, RPROJ, X(3*NATOMS), DUMMYGRAD(3*NATOMS), DUMMYE, THETA2, FRPISQ, OPOTEL
      PARAMETER (PI=3.141592654D0)
      INTEGER J1, J2, JMAX, NP, J3, JMAX2, RANATOM, INDEXD(NATOMS), NCOOPDONE, ISTAT, NDUMMY, NTRIES, NTRIESMAX



      IF (PERMOPT.AND.PERIODIC) RETURN 
      
      IF (MODEL1T) THEN
         RANDOM=(DPRAND()-0.5D0)*2.0D0
         RETURN
      ENDIF
      FRPISQ = 4.D0*PI*PI
      NTRIESMAX=100

      IF ((.NOT.NORESET).AND.(.NOT.PERMOPT).AND.(.NOT.DIFFRACTT).AND.(.NOT.BLNT).AND.(.NOT.PERIODIC)
     &     .AND.(.NOT.GAUSST).AND.(.NOT.(CSMT.AND.(.NOT.SYMMETRIZECSM))).AND.(.NOT.PERCOLATET)) THEN
         DO J1=1,NATOMS
            IF ((.NOT.RIGID).OR.(J1.LE.NATOMS/2)) THEN
               J2=3*J1
               DUMMY2=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
               IF (DUMMY2.GT.RADIUS) THEN
                  IF (AMBERT) THEN ! jmc49 We don't really want a container at all in amber9, but this bit of code is being used 
                     WRITE(MYUNIT,'(A,I5,5F20.10)')'J1,RAD,D2,x,y,z=',J1,RADIUS,DUMMY2,COORDS(J2-2,NP),COORDS(J2-1,NP),COORDS(J2,NP)
                     WRITE(MYUNIT,'(A,I6)') 'WARNING initial coordinate outside container for particle ',J1
                     WRITE(MYUNIT,'(A)') 'NOT repositioning particle for AMBER9'
                  END IF
                  WRITE(MYUNIT,'(A,I5,5F20.10)') 'J1,RAD,D2,x,y,z=',J1,RADIUS,DUMMY2,COORDS(J2-2,NP),COORDS(J2-1,NP),COORDS(J2,NP)
                  WRITE(MYUNIT,'(A,I6)') 'WARNING initial coordinate outside container - reposition particle ',J1
                  IF (DEBUG) THEN
                     WRITE(77,*) NATOMS
                     WRITE(77,*) ' '
                     WRITE(77,'(A3,3F20.10)') ('LA ',COORDS(3*(J3-1)+1,NP),COORDS(3*(J3-1)+2,NP),COORDS(3*(J3-1)+3,NP),J3=1,NATOMS)
                  ENDIF

               ENDIF
            ENDIF
         ENDDO
      ENDIF

      IF (WENZEL) THEN
11       RANDOM=(DPRAND()-0.5D0)*2.0D0
         IF ((COORDS(1,NP).GT.1.0D0).OR.(COORDS(1,NP).LT.0.0D0)) GOTO 11
12       RANDOM=(DPRAND()-0.5D0)*2.0D0
         IF ((COORDS(2,NP).GT.1.0D0).OR.(COORDS(2,NP).LT.0.0D0)) GOTO 12
         RETURN
      ELSE IF (PERMOPT.OR.(CSMT.AND.(.NOT.SYMMETRIZECSM))) THEN
         RANDOM=(DPRAND()-0.5D0)*2.0D0
         ANGLE=RANDOM*STEP(NP) ! in radians
         SINT=SIN(ANGLE)
         DO J1=1,NATOMS
            TY= COST*COORDS(3*(J1-1)+2,NP)+SINT*COORDS(3*(J1-1)+3,NP)
            TZ=-SINT*COORDS(3*(J1-1)+2,NP)+COST*COORDS(3*(J1-1)+3,NP)
         ENDDO
         RANDOM=(DPRAND()-0.5D0)*2.0D0
         ANGLE=RANDOM*STEP(NP) ! in radians
         SINT=SIN(ANGLE)
         DO J1=1,NATOMS
            TX= COST*COORDS(3*(J1-1)+1,NP)+SINT*COORDS(3*(J1-1)+3,NP)
            TZ=-SINT*COORDS(3*(J1-1)+1,NP)+COST*COORDS(3*(J1-1)+3,NP)
         ENDDO
         RANDOM=(DPRAND()-0.5D0)*2.0D0
         ANGLE=RANDOM*STEP(NP) ! in radians
         SINT=SIN(ANGLE)
         DO J1=1,NATOMS
            TX= COST*COORDS(3*(J1-1)+1,NP)+SINT*COORDS(3*(J1-1)+2,NP)
            TY=-SINT*COORDS(3*(J1-1)+1,NP)+COST*COORDS(3*(J1-1)+2,NP)
         ENDDO

         RETURN
      ENDIF
      XMASS=0.0D0; YMASS=0.0D0; ZMASS=0.0D0
      DO J1=1,NATOMS
         XMASS=XMASS+COORDS(3*(J1-1)+1,NP)
         YMASS=YMASS+COORDS(3*(J1-1)+2,NP)
         ZMASS=ZMASS+COORDS(3*(J1-1)+3,NP)
      ENDDO
      XMASS=XMASS/NATOMS; YMASS=YMASS/NATOMS; ZMASS=ZMASS/NATOMS

      DMAX=-1.0D0
      VMAX=-1.0D6
      VMAX2=-1.0D6
      VMIN=1.0D6
      DO J1=1,NATOMS
         J2=3*J1
         DIST(J1)= DSQRT( COORDS(J2-2,NP)**2+        COORDS(J2-1,NP)**2+        COORDS(J2,NP)**2)
         IF ((CMDIST(J1).GT.CMMAX).AND.(J1.LE.NATOMS-NCORE(NP))) CMMAX=CMDIST(J1)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1)
            IF (VAT(J1,NP).GT.VMAX) THEN
               VMAX=VAT(J1,NP)
               JMAX=J1
            ELSE IF ((VAT(J1,NP).LT.VMAX).AND.(VAT(J1,NP).GT.VMAX2)) THEN
               VMAX2=VAT(J1,NP)
               JMAX2=J1
            ENDIF
         IF (VAT(J1,NP).LT.VMIN) VMIN=VAT(J1,NP)
      ENDDO
  
      IF (DECAY) THEN
9        RANATOM=NINT(0.5D0+NATOMS*DPRAND())
         IF (RANATOM.EQ.JMAX) GOTO 9 ! don't choose the atom that might undergo a surface move
         WRITE(MYUNIT,'(A,I6)') 'atom undergoing maximum displacement is number ',RANATOM
         DO J1=1,NATOMS-NSEED
            DUMMY=((COORDS(3*J1-2,NP)-COORDS(3*RANATOM-2,NP))**2+
     1             (COORDS(3*J1-1,NP)-COORDS(3*RANATOM-1,NP))**2+
     2             (COORDS(3*J1,NP)-  COORDS(3*RANATOM,NP))**2)
            EXPDECAY(J1)=EXP(-DECAYPARAM*DUMMY)
         ENDDO
      ENDIF
      IF (NCORE(NP).GT.0) THEN
         NSURFMOVES(NP)=NSURFMOVES(NP)+1
         XC=0.0D0; YC=0.0D0; ZC=0.0D0
         DO J1=NATOMS,NATOMS-NCORE(NP)+1,-1
            XC=XC+COORDS(3*(J1-1)+1,NP)
            YC=YC+COORDS(3*(J1-1)+2,NP)
            ZC=ZC+COORDS(3*(J1-1)+3,NP)
         ENDDO
         XC=XC/NCORE(NP); YC=YC/NCORE(NP); ZC=ZC/NCORE(NP)
         IF (DEBUG) WRITE(MYUNIT,'(A,3F12.4)') 'takestep> centre of coordinates for frozen atoms: ',XC, YC, ZC
         IF (MOVESHELLT.AND.(DPRAND().GT.(1.0D0-SHELLPROB))) THEN 
            WRITE(MYUNIT,'(A,I1,A,I8)') '[',NP,']takestep> shell move number ',NSURFMOVES(NP)

            RVEC(1)=(DPRAND()-0.5D0)*2.0D0
            RVEC(2)=(DPRAND()-0.5D0)*2.0D0
            RVEC(3)=(DPRAND()-0.5D0)*2.0D0
            DUMMY=SQRT(RVEC(1)**2+RVEC(2)**2+RVEC(3)**2)
            RVEC(1)=RVEC(1)/DUMMY; RVEC(2)=RVEC(2)/DUMMY; RVEC(3)=RVEC(3)/DUMMY
            ANGLE=DPRAND()*PI*2.0D0
            SINT=SIN(ANGLE)

            IF (DEBUG) WRITE(MYUNIT,'(A,I1,A,F10.2,A,3F12.4)') '[',NP,']takestep> angle=',ANGLE,' axis: ',RVEC(1:3)
            DO J1=1,NATOMS-NCORE(NP)
               XL=COORDS(3*(J1-1)+1,NP); YL=COORDS(3*(J1-1)+2,NP); ZL=COORDS(3*(J1-1)+3,NP)
               DUMMY=SQRT((XL-XC)**2+(YL-YC)**2+(ZL-ZC)**2)
   
               RDOTN=(XL-XC)*RVEC(1)+(YL-YC)*RVEC(2)+(ZL-ZC)*RVEC(3)
               TX=(XL-XC)*COST + RVEC(1)*RDOTN*(1.0D0-COST)-((YL-YC)*RVEC(3)-(ZL-ZC)*RVEC(2))*SINT
               TY=(YL-YC)*COST + RVEC(2)*RDOTN*(1.0D0-COST)-((ZL-ZC)*RVEC(1)-(XL-XC)*RVEC(3))*SINT
               TZ=(ZL-ZC)*COST + RVEC(3)*RDOTN*(1.0D0-COST)-((XL-XC)*RVEC(2)-(YL-YC)*RVEC(1))*SINT
               IF (DUMMY.GT.0.1D0) THEN
               ELSE 
               ENDIF
            ENDDO
            IF (NSURFMOVES(NP).GE.SHELLMOVEMAX) THEN
               SHELLMOVES(NP)=.FALSE.
               NCORE(NP)=0
            ENDIF
            RETURN
         ENDIF
333      CONTINUE
      ENDIF

      IF (COOP.AND.(NCORE(NP).GE.MAX(2,NCOOP))) THEN
8        IF (NCORE(NP).GT.0) THEN
            RANATOM=NATOMS-NCORE(NP)+NINT(0.5D0+NCORE(NP)*DPRAND())
         ELSE
            RANATOM=NINT(0.5D0+(NATOMS-NSEED)*DPRAND())
         ENDIF
         IF (RANATOM.EQ.JMAX) GOTO 8 ! don't choose the atom that might undergo a surface move
         IF (DEBUG) WRITE(MYUNIT,'(A,I6)') 'takestep> randomly selected atom for coop move is number ',RANATOM
         DO J1=1,NATOMS-NSEED

     1                 (COORDS(3*J1-1,NP)-COORDS(3*RANATOM-1,NP))**2+
     2                 (COORDS(3*J1,NP)-  COORDS(3*RANATOM,NP))**2)
            IF ((NCORE(NP).GT.0).AND.(J1.GT.NATOMS-NCORE(NP))) CDIST(J1)=1.0D100
            INDEXD(J1)=J1
         ENDDO
         RANDOMX=(DPRAND()-0.5D0)*2.0D0
         RANDOMY=(DPRAND()-0.5D0)*2.0D0
         RANDOMZ=(DPRAND()-0.5D0)*2.0D0
         DUMMY2=SQRT(RANDOMX**2+RANDOMY**2+RANDOMZ**2)
         RANDOMX=RANDOMX/DUMMY2
         RANDOMY=RANDOMY/DUMMY2
         RANDOMZ=RANDOMZ/DUMMY2
         NCOOPDONE=0
      ENDIF
      DO J1=1,NATOMS-NSEED
         IF (FROZEN(J1)) CYCLE
         IF ((COREFRAC.EQ.0.0D0).AND.(J1.GT.NATOMS-NCORE(NP))) CYCLE ! no point taking a zero step
         NTRIES=0
10       J2=3*J1
         IF ((NCORE(NP).GT.0).AND.(J1.GT.NATOMS-NCORE(NP))) THEN
            LOCALSTEP=STEP(NP)*COREFRAC ! smaller step for core
         ELSE
            LOCALSTEP=STEP(NP)
         ENDIF
         IF (RIGID.AND.(J1.GT.NATOMS/2)) THEN
            LOCALSTEP=0.0D0
            IF (OMOVE(NP)) LOCALSTEP=OSTEP(NP)
         ELSE IF (RIGID.AND.(J1.LE.NATOMS/2)) THEN
            LOCALSTEP=0.0D0
            IF (TMOVE(NP)) LOCALSTEP=STEP(NP)
         ENDIF
         IF (((VAT(J1,NP).GT.ASTEP(NP)*VMIN).AND.(J1.EQ.JMAX)).AND.(.NOT.RIGID).AND.(.NOT.BLNT).AND. 
     &         (.NOT.DIFFRACTT).AND.(.NOT.GAUSST).AND.(.NOT.PERCOLATET) 
     &        .AND.(.NOT.NORESET).AND.(.NOT.PERIODIC).AND.(.NOT.THOMSONT)
     &        .AND.(.NOT.((NCORE(NP).GT.0).AND.(J1.GT.NATOMS-NCORE(NP))))) THEN

            IF (DEBUG) WRITE(MYUNIT,'(A,I4,A,F12.4,A,F12.4,A,I4,A,F12.4)') 'angular move for atom ',J1, 
     &           ' V=',VMAX,' Vmin=',VMIN,' next most weakly bound atom is ',JMAX2,' V=',VMAX2

           THETA=DPRAND()*PI
           PHI=DPRAND()*PI*2.0D0
           DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
           IF (DUMMY.GT.RADIUS) THEN
              DUMMY=SQRT(RADIUS*0.99D0/DUMMY)
           ENDIF
         ELSE IF ((NATOMS-NSEED.EQ.1).AND.(NATOMS.GT.1)) THEN 
           IF (DEBUG) WRITE(MYUNIT,'(A,I4,A,F12.4,A,2F12.4)') 
     1                'angular move for atom ',J1,' V=',VAT(J1,NP),' Vmin, Vmax=',VMIN,VMAX
           THETA=DPRAND()*PI
           PHI=DPRAND()*PI*2.0D0
           DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
           IF (DUMMY.GT.RADIUS) THEN
              DUMMY=SQRT(RADIUS*0.99D0/DUMMY)
           ENDIF
         ELSE IF (DECAY) THEN
            RANDOMX=(DPRAND()-0.5D0)*2.0D0
            RANDOMY=(DPRAND()-0.5D0)*2.0D0
            RANDOMZ=(DPRAND()-0.5D0)*2.0D0
            DUMMY2=SQRT(RANDOMX**2+RANDOMY**2+RANDOMZ**2)
         ELSE IF (COOP.AND.(NCORE(NP).GE.MAX(2,NCOOP))) THEN
            J2=3*INDEXD(J1)
            DUMMY=CDIST(J1)
            IF (NCORE(NP).GT.0) THEN
               J2=3*INDEXD(J1)
               DUMMY=CDIST(J1)
            ENDIF
            IF ((NCOOPDONE.LE.NCOOP+1).AND.(DUMMY.LT.COOPCUT)) THEN
               NCOOPDONE=NCOOPDONE+1
               WRITE(MYUNIT,'(A,I6,A,F12.4)') 'takestep> coop move for atom ',J2/3,' cdist=',DUMMY
            ELSE
            ENDIF
         ELSE IF (.NOT.FIXD) THEN
           RANDOM=DPRAND()
           IF ((VMIN-VMAX.EQ.0.0D0).OR.(EFAC.EQ.0.0D0)) THEN
                 IF (NCORE(NP).GT.0) THEN ! project out radial component and rescale to same length
                    DELX=(DPRAND()-0.5D0)*2.0D0*LOCALSTEP
                    DELY=(DPRAND()-0.5D0)*2.0D0*LOCALSTEP
                    DELZ=(DPRAND()-0.5D0)*2.0D0*LOCALSTEP
                    IF (J1.LE.NATOMS-NCORE(NP)) THEN
                       SLENGTH=SQRT(DELX**2+DELY**2+DELZ**2)
                       DUMMY=SQRT((COORDS(J2-2,NP)-XC)**2+(COORDS(J2-1,NP)-YC)**2+(COORDS(J2,NP)-ZC)**2)
                       RPROJ=(COORDS(J2-2,NP)-XC)*DELX+(COORDS(J2-1,NP)-YC)*DELY+(COORDS(J2,NP)-ZC)*DELZ
                       IF (DUMMY.NE.0.0D0) THEN
                          RPROJ=RPROJ/DUMMY
                          DELX=DELX-RPROJ*(COORDS(J2-2,NP)-XC)/DUMMY
                          DELY=DELY-RPROJ*(COORDS(J2-1,NP)-YC)/DUMMY
                          DELZ=DELZ-RPROJ*(COORDS(J2,NP)-ZC)/DUMMY
                          DUMMY=SQRT(DELX**2+DELY**2+DELZ**2)
                          IF (DUMMY.NE.0.0D0) THEN
                             DELX=DELX*SLENGTH/DUMMY
                             DELY=DELY*SLENGTH/DUMMY
                             DELZ=DELZ*SLENGTH/DUMMY
                          ENDIF
                       ENDIF
                    ENDIF
                 ELSE
                 RANDOM=(DPRAND()-0.5D0)*2.0D0
                 IF (GAUSST) RANDOM=RANDOM/GKSMALL(J2-2) ! scale gauss steps by 1/GKSMALL
                 RANDOM=(DPRAND()-0.5D0)*2.0D0
                 IF (GAUSST) RANDOM=RANDOM/GKSMALL(J2-1)
                 RANDOM=(DPRAND()-0.5D0)*2.0D0
                 IF (GAUSST) RANDOM=RANDOM/GKSMALL(J2)
              ENDIF
           ELSE 
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              DUMMY=1.0D0+EAMP*TANH(-2.0D0*EFAC*(VAT(J1,NP)-(VMAX+VMIN)/2.0D0)/(VMAX-VMIN))
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              RANDOM=(DPRAND()-0.5D0)*2.0D0
           ENDIF
           IF ((.NOT.PERIODIC).AND.(.NOT.AMBERT).AND.(.NOT.(RIGID.AND.((J1.GT.NATOMS/2)))).AND.(.NOT.BLNT)
     1     .AND.(.NOT.PERCOLATET).AND.(.NOT.DIFFRACTT).AND.(.NOT.THOMSONT).AND.(.NOT.GAUSST)) THEN
              DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
              IF (DUMMY.GT.RADIUS) THEN
                 NTRIES=NTRIES+1
                 IF (NTRIES.GT.NTRIESMAX) THEN
                    WRITE(MYUNIT,'(A,I6,A)') 'takestep> WARNING atom ',J1,' persistently outside container'
                 ELSE
                    IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,I6)') 'takestep> WARNING atom ',J1,' outside container, NTRIES=',NTRIES
                    GOTO 10
                 ENDIF
              ENDIF
           ENDIF
         ENDIF
         IF (TOSI.OR.WELCH) THEN
            DO J3=1,J1-1
               DUMMY=(COORDS(J2-2,NP)-COORDS(3*(J3-1)+1,NP))**2
     1              +(COORDS(J2-1,NP)-COORDS(3*(J3-1)+2,NP))**2
     2              +(COORDS(J2  ,NP)-COORDS(3*(J3-1)+3,NP))**2
               IF (DUMMY.LT.1.0D0) GOTO 10 
            ENDDO
         ENDIF
      ENDDO

      IF (FIXD) CALL HSMOVE(COORDS(1:3*NATOMS,1:NPAR),NP,NHSMOVE)
      IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,NP))
      IF (FIXCOM.AND.(.NOT.SEEDT)) CALL CENTRECOM(COORDS(:3*NATOMS,NP))

      IF (NSYMREM.GT.0) THEN
      ENDIF

      IF (DBPT.OR.DBPTDT.OR.DMBLMT.OR.LWOTPT.OR.MSSTOCKT.OR.MSTBINT.OR.NCAPT.OR.NPAHT.OR.NTIPT.OR.MULTPAHAT.OR.PAHAT.OR.STOCKAAT
     &    .OR.PAHW99T.OR.TDHDT.OR.PATCHY.OR.PAPT.OR.SILANET) THEN

         DO J1 = NATOMS/2 + 1, NATOMS
            J2      = 3*J1
            THETA2  = DOT_PRODUCT(COORDS(J2-2:J2,NP),COORDS(J2-2:J2,NP))
            IF (THETA2 > FRPISQ) THEN
               THETA2   = DSQRT(THETA2)
               THETA    = THETA2 - INT(THETA2/(2.D0*PI))*2.D0*PI
            ENDIF
         ENDDO

      ENDIF 

      RETURN
      END
      SUBROUTINE KEEPSYM(NP)
      USE commons
      IMPLICIT NONE
      INTEGER NP, J2, SYMINDEX(NATOMS), J3, J4
      DOUBLE PRECISION LCOORDS(3*NATOMS), NEWQ(3*NATOMS), SYMDELTA(3*NATOMS), DELTA(3*NATOMS),  SYMOP1(3,3)
      DOUBLE PRECISION STEPLENGTH, NEWSTEPLENGTH, ODIST1, ODIST2, DMIN, DUMMY
      LOGICAL ASSIGNED(NATOMS), BAD

      LCOORDS(1:3*NATOMS)=COORDSO(1:3*NATOMS,NP)
      DELTA(1:3*NATOMS)=COORDS(1:3*NATOMS,NP)-COORDSO(1:3*NATOMS,NP)
      SYMDELTA(1:3*NATOMS)=DELTA(1:3*NATOMS)

      DO J2=1,NSYMREM
         BAD=.FALSE.
         SYMOP1(1:3,1:3)=SYMREM(J2,1:3,1:3)
         ASSIGNED(1:NATOMS)=.FALSE.
         SYMINDEX(1:NATOMS)=0
         DO J3=1,NATOMS
            DMIN=1.0D100
            DO J4=1,NATOMS
               DUMMY=(LCOORDS(3*(J4-1)+1)-NEWQ(3*(J3-1)+1))**2+
     &               (LCOORDS(3*(J4-1)+2)-NEWQ(3*(J3-1)+2))**2+
     &               (LCOORDS(3*(J4-1)+3)-NEWQ(3*(J3-1)+3))**2
               IF (DUMMY.LT.DMIN) THEN
                  IF (ASSIGNED(J4)) THEN
                  ELSE
                     IF (SYMINDEX(J3).GT.0) ASSIGNED(SYMINDEX(J3))=.FALSE.
                     SYMINDEX(J3)=J4
                     ASSIGNED(J4)=.TRUE.
                     DMIN=DUMMY
                  ENDIF
               ENDIF
            ENDDO
            IF (DEBUG.AND.(SQRT(DMIN).GT.SYMTOL5)) THEN
               WRITE(MYUNIT, '(2(A,I5),A,F12.3)') 'WARNING closest image to atom ',J3,' is ',SYMINDEX(J3), 
     &                                         ' distance=',SQRT(DMIN)
               BAD=.TRUE.
            ENDIF
         ENDDO

         DO J3=1,NATOMS
            SYMDELTA(3*(J3-1)+1)=SYMDELTA(3*(J3-1)+1)+NEWQ(3*(SYMINDEX(J3)-1)+1)
            SYMDELTA(3*(J3-1)+2)=SYMDELTA(3*(J3-1)+2)+NEWQ(3*(SYMINDEX(J3)-1)+2)
            SYMDELTA(3*(J3-1)+3)=SYMDELTA(3*(J3-1)+3)+NEWQ(3*(SYMINDEX(J3)-1)+3)
         ENDDO
      ENDDO

      SYMDELTA(1:3*NATOMS)=SYMDELTA(1:3*NATOMS)/(1+NSYMREM)

      DO J2=1,NATOMS
         ODIST1=SUM(COORDSO(3*(J2-1)+1:3*(J2-1)+3,NP)**2)
         ODIST2=SUM((COORDSO(3*(J2-1)+1:3*(J2-1)+3,NP)+SYMDELTA(3*(J2-1)+1:3*(J2-1)+3))**2)
         IF (ODIST2.NE.0.0D0) COORDS(3*(J2-1)+1:3*(J2-1)+3,NP)=(COORDSO(3*(J2-1)+1:3*(J2-1)+3,NP)
     &        +SYMDELTA(3*(J2-1)+1:3*(J2-1)+3)) * SQRT(ODIST1/ODIST2)
      ENDDO


    
      RETURN
      END
