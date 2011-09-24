
      SUBROUTINE TAKESTEP(NP)

      USE COMMONS
      USE V
      USE F

      IMPLICIT NONE
      
      INTEGER,INTENT(IN) ::    NP 
      DOUBLE PRECISION ::   RANDOM
      DOUBLE PRECISION ::   DMAX,VMAX,VMIN

!      DOUBLE PRECISION RANDOM, XMASS, YMASS, ZMASS, LOCALSTEP, DUMMY2, CDIST(NATOMS), RDOTN, XL, YL, ZL,&
     !&                 DIST(3*NATOMS), DMAX, VMAX, VMIN, VMAX2, EXPDECAY(NATOMS), XC, YC, ZC, ANGLE, COST, SINT,&
     !&                 THETA, PHI, PI, DUMMY, CMDIST(NATOMS), CMMAX, RANDOMX, RANDOMY, RANDOMZ, RVEC(3), TX, TY, TZ,&
     !&                 DELX, DELY, DELZ, SLENGTH, RPROJ, X(3*NATOMS), DUMMYGRAD(3*NATOMS), DUMMYE, THETA2, FRPISQ, OPOTEL
      !PARAMETER (PI=3.141592654D0)
      !INTEGER J1, J2, JMAX, J3, JMAX2, RANATOM, INDEXD(NATOMS), NCOOPDONE, ISTAT, NDUMMY, NTRIES, NTRIESMAX
!!     INTEGER NQTOT
!     COMMON /TOT/ NQTOT

!  Calling CENTRE if NORESET is .TRUE. can lead to problems with COORDSO containing an atom
!  outside the permitted radius. Then it may be impossible to take a step that keeps all the
!  atoms inside.
!
      !FRPISQ = 4.D0*PI*PI
      !NTRIESMAX=100

! COORDSO and VATO should now contain the correct COORDS and VAT even if
! these entries are changed by SYMMETRIZE or NEWRESET.
!
       IF (ABS(COORDSO(1,NP)-COORDS(1,NP)).GT.1.0D-3) THEN
          WRITE(LFH,'(A,2G20.10)'),'takestep> WARNING - coordso will be changed: ',COORDSO(1,NP),COORDS(1,NP)
       ENDIF
       DO J1=1,3*(NATOMS-NSEED)
          COORDSO(J1,NP)=COORDS(J1,NP)
       ENDDO
       DO J1=1,NATOMS
          VATO(J1,NP)=VAT(J1,NP)
       ENDDO

!	csw34> FIND CENTRE OF COORDINATES (should be made a function!)

      XMASS=0.0D0; YMASS=0.0D0; ZMASS=0.0D0
      DO J1=1,NATOMS
         XMASS=XMASS+COORDS(3*(J1-1)+1,NP)
         YMASS=YMASS+COORDS(3*(J1-1)+2,NP)
         ZMASS=ZMASS+COORDS(3*(J1-1)+3,NP)
      ENDDO
      XMASS=XMASS/NATOMS; YMASS=YMASS/NATOMS; ZMASS=ZMASS/NATOMS
!
!  Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
!  and the pair energy of the most tightly bound atom, VMIN. An angular step is
!  taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
!  DMAX (or CMMAX from CM of the cluster).
!

      DMAX=-1.0D0
      VMAX=-1.0D6
      VMAX2=-1.0D6
      VMIN=1.0D6
      CMMAX=-1.0D0
      DO J1=1,NATOMS
         J2=3*J1
         DIST(J1)= DSQRT( COORDS(J2-2,NP)**2+        COORDS(J2-1,NP)**2+        COORDS(J2,NP)**2)
         CMDIST(J1)=SQRT((COORDS(J2-2,NP)-XMASS)**2+(COORDS(J2-1,NP)-YMASS)**2+(COORDS(J2,NP)-ZMASS)**2)
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
!!{{{
!  If DECAY is true then select an atom at random, move this one randomly
!  by the maximum allowed amount, and move the others in random directions
!  of decaying magnitude, depending on how far they are from the chosen atom.
!
  
      !IF (DECAY) THEN
!9        RANATOM=NINT(0.5D0+NATOMS*DPRAND())
         !IF (RANATOM.EQ.JMAX) GOTO 9 ! don't choose the atom that might undergo a surface move
         !WRITE(LFH,'(A,I6)') 'atom undergoing maximum displacement is number ',RANATOM
         !DO J1=1,NATOMS-NSEED
            !DUMMY=((COORDS(3*J1-2,NP)-COORDS(3*RANATOM-2,NP))**2+ &
     !&             (COORDS(3*J1-1,NP)-COORDS(3*RANATOM-1,NP))**2+ &
     !&             (COORDS(3*J1,NP)-  COORDS(3*RANATOM,NP))**2) 
            !EXPDECAY(J1)=EXP(-DECAYPARAM*DUMMY)
         !ENDDO
      !ENDIF
!
!  If MOVESHELL is true then we try a random angular move for all the frozen atoms 
!  about the centre of coordinates of the frozen set with probability SHELLPROB.
!
!     IF (SHELLMOVES(NP)) THEN
      !IF (NCORE(NP).GT.0) THEN
         !NSURFMOVES(NP)=NSURFMOVES(NP)+1
         !XC=0.0D0; YC=0.0D0; ZC=0.0D0
         !DO J1=NATOMS,NATOMS-NCORE(NP)+1,-1
            !XC=XC+COORDS(3*(J1-1)+1,NP)
            !YC=YC+COORDS(3*(J1-1)+2,NP)
            !ZC=ZC+COORDS(3*(J1-1)+3,NP)
         !ENDDO
         !XC=XC/NCORE(NP); YC=YC/NCORE(NP); ZC=ZC/NCORE(NP)
         !IF (DEBUG) WRITE(LFH,'(A,3F12.4)') 'takestep> centre of coordinates for frozen atoms: ',XC, YC, ZC
         !IF (MOVESHELLT.AND.(DPRAND().GT.(1.0D0-SHELLPROB))) THEN 
            !WRITE(LFH,'(A,I1,A,I8)') '[',NP,']takestep> shell move number ',NSURFMOVES(NP)

            !RVEC(1)=(DPRAND()-0.5D0)*2.0D0
            !RVEC(2)=(DPRAND()-0.5D0)*2.0D0
            !RVEC(3)=(DPRAND()-0.5D0)*2.0D0
            !DUMMY=SQRT(RVEC(1)**2+RVEC(2)**2+RVEC(3)**2)
            !RVEC(1)=RVEC(1)/DUMMY; RVEC(2)=RVEC(2)/DUMMY; RVEC(3)=RVEC(3)/DUMMY
            !ANGLE=DPRAND()*PI*2.0D0
            !COST=COS(ANGLE)
            !SINT=SIN(ANGLE)

            !IF (DEBUG) WRITE(LFH,'(A,I1,A,F10.2,A,3F12.4)') '[',NP,']takestep> angle=',ANGLE,' axis: ',RVEC(1:3)
!!
!!  Rotate all the non-core atoms through ANGLE about RVEC. Use rotation formula
!!  from Goldstein p. 165.
!!  
!!           DO J1=NCORE(NP)+1,NATOMS
            !DO J1=1,NATOMS-NCORE(NP)
               !XL=COORDS(3*(J1-1)+1,NP); YL=COORDS(3*(J1-1)+2,NP); ZL=COORDS(3*(J1-1)+3,NP)
               !DUMMY=SQRT((XL-XC)**2+(YL-YC)**2+(ZL-ZC)**2)
   
               !RDOTN=(XL-XC)*RVEC(1)+(YL-YC)*RVEC(2)+(ZL-ZC)*RVEC(3)
               !TX=(XL-XC)*COST + RVEC(1)*RDOTN*(1.0D0-COST)-((YL-YC)*RVEC(3)-(ZL-ZC)*RVEC(2))*SINT
               !TY=(YL-YC)*COST + RVEC(2)*RDOTN*(1.0D0-COST)-((ZL-ZC)*RVEC(1)-(XL-XC)*RVEC(3))*SINT
               !TZ=(ZL-ZC)*COST + RVEC(3)*RDOTN*(1.0D0-COST)-((XL-XC)*RVEC(2)-(YL-YC)*RVEC(1))*SINT
               !IF (DUMMY.GT.0.1D0) THEN
                  !COORDS(3*(J1-1)+1,NP)=(XC+TX)*(DUMMY+1.0D0)/DUMMY
                  !COORDS(3*(J1-1)+2,NP)=(YC+TY)*(DUMMY+1.0D0)/DUMMY
                  !COORDS(3*(J1-1)+3,NP)=(ZC+TZ)*(DUMMY+1.0D0)/DUMMY
               !ELSE 
                  !COORDS(3*(J1-1)+1,NP)=(XC+TX)
                  !COORDS(3*(J1-1)+2,NP)=(YC+TY)
                  !COORDS(3*(J1-1)+3,NP)=(ZC+TZ)
               !ENDIF
            !ENDDO
            !IF (NSURFMOVES(NP).GE.SHELLMOVEMAX) THEN
               !SHELLMOVES(NP)=.FALSE.
               !NCORE(NP)=0
            !ENDIF
            !RETURN
         !ENDIF
!!        IF (NSURFMOVES(NP).GE.SHELLMOVEMAX) THEN
!!           NSURFMOVES(NP)=-1
!!           SHELLMOVES(NP)=.FALSE.
!!           NCORE(NP)=0
!!        ENDIF
!333      CONTINUE
      !ENDIF

!     IF (NSURFMOVES(NP).LT.0) THEN
!        NSURFMOVES(NP)=NSURFMOVES(NP)-1
!        IF (NSURFMOVES(NP).LE.-SHELLMOVEMAX) THEN
!           NSURFMOVES(NP)=0 
!        ENDIF
!     ENDIF
!
!  If COOP is true then select an atom at random, move this one randomly
!  by the maximum allowed amount, and move its NCOOP nearest neighbours by
!  the same displacement.
!  Note that the core atoms are ordered LAST, not FIRST !!
!
!     IF (COOP.AND.(NATOMS-NCORE(NP).GE.MAX(2,NCOOP))) THEN
      !IF (COOP.AND.(NCORE(NP).GE.MAX(2,NCOOP))) THEN
!!     IF (COOP) THEN
!8        IF (NCORE(NP).GT.0) THEN
!!           RANATOM=NINT(0.5D0+(NATOMS-NCORE(NP))*DPRAND())
!!           RANATOM=NINT(0.5D0+(NATOMS-NSEED)*DPRAND())
            !RANATOM=NATOMS-NCORE(NP)+NINT(0.5D0+NCORE(NP)*DPRAND())
         !ELSE
            !RANATOM=NINT(0.5D0+(NATOMS-NSEED)*DPRAND())
         !ENDIF
         !IF (RANATOM.EQ.JMAX) GOTO 8 ! don't choose the atom that might undergo a surface move
         !IF (DEBUG) WRITE(LFH,'(A,I6)') 'takestep> randomly selected atom for coop move is number ',RANATOM
         !DO J1=1,NATOMS-NSEED

            !CDIST(J1)=((COORDS(3*J1-2,NP)-COORDS(3*RANATOM-2,NP))**2+&
     !&                 (COORDS(3*J1-1,NP)-COORDS(3*RANATOM-1,NP))**2+&
     !&                 (COORDS(3*J1,NP)-  COORDS(3*RANATOM,NP))**2)
!!           IF (SHELLMOVES(NP).AND.(J1.LE.NCORE(NP))) CDIST(J1)=1.0D100
!!           IF (SHELLMOVES(NP).AND.(J1.GT.NATOMS-NCORE(NP))) CDIST(J1)=1.0D100
            !IF ((NCORE(NP).GT.0).AND.(J1.GT.NATOMS-NCORE(NP))) CDIST(J1)=1.0D100
            !INDEXD(J1)=J1
         !ENDDO
         !CALL SORT4(NATOMS-NSEED,NATOMS,CDIST,INDEXD)
!!        PRINT *,'INDEXD=',INDEXD(1:NATOMS)
!!        PRINT *,'CDIST=',CDIST(1:NATOMS)
!!        CALL FLUSH(6,ISTAT)
         !RANDOMX=(DPRAND()-0.5D0)*2.0D0
         !RANDOMY=(DPRAND()-0.5D0)*2.0D0
         !RANDOMZ=(DPRAND()-0.5D0)*2.0D0
         !DUMMY2=SQRT(RANDOMX**2+RANDOMY**2+RANDOMZ**2)
         !RANDOMX=RANDOMX/DUMMY2
         !RANDOMY=RANDOMY/DUMMY2
         !RANDOMZ=RANDOMZ/DUMMY2
         !NCOOPDONE=0
      !ENDIF
!!
!	csw34> For the below loop, J1 is the atom number, and J2 is set to point to 3*J1, meaning
!	the z-coordinate of that atom. The y coordinate is then J2-1 and x is J2-2. 
!!}}}

      DO J1=1,NATOMS-NSEED
         NTRIES=0
10       J2=3*J1
         LOCALSTEP=STEP(NP)

         !IF ((NCORE(NP).GT.0).AND.(J1.GT.NATOMS-NCORE(NP))) THEN
            !LOCALSTEP=STEP(NP)*COREFRAC ! smaller step for core
         !ELSE
            !LOCALSTEP=STEP(NP)
         !ENDIF
        ! IF (RIGID.AND.(J1.GT.NATOMS/2)) THEN
            !LOCALSTEP=0.0D0
            !IF (OMOVE(NP)) LOCALSTEP=OSTEP(NP)
         !ELSE IF (RIGID.AND.(J1.LE.NATOMS/2)) THEN
            !LOCALSTEP=0.0D0
            !IF (TMOVE(NP)) LOCALSTEP=STEP(NP)
         !ENDIF
         !ELSE IF ((NATOMS-NSEED.EQ.1).AND.(NATOMS.GT.1)) THEN 
!         IF ((NATOMS-NSEED.EQ.1).AND.(NATOMS.GT.1)) THEN 
           !IF (DEBUG) WRITE(LFH,'(A,I4,A,F12.4,A,2F12.4)') &
     !&                'angular move for atom ',J1,' V=',VAT(J1,NP),' Vmin, Vmax=',VMIN,VMAX
           !THETA=DPRAND()*PI
           !PHI=DPRAND()*PI*2.0D0
           !COORDS(J2-2,NP)=XMASS+(CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
           !COORDS(J2-1,NP)=YMASS+(CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
           !COORDS(J2,NP)=  ZMASS+(CMMAX+1.0D0)*DCOS(THETA)
           !DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
           !IF (DUMMY.GT.RADIUS) THEN
              !DUMMY=SQRT(RADIUS*0.99D0/DUMMY)
              !COORDS(J2-2,NP)=COORDS(J2-2,NP)*DUMMY
              !COORDS(J2-1,NP)=COORDS(J2-1,NP)*DUMMY
              !COORDS(J2,NP)=COORDS(J2,NP)*DUMMY
           !ENDIF
         !ELSE IF (DECAY) THEN
            !RANDOMX=(DPRAND()-0.5D0)*2.0D0
            !RANDOMY=(DPRAND()-0.5D0)*2.0D0
            !RANDOMZ=(DPRAND()-0.5D0)*2.0D0
            !DUMMY2=SQRT(RANDOMX**2+RANDOMY**2+RANDOMZ**2)
            !COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOMX*EXPDECAY(J1)/DUMMY2
            !COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOMY*EXPDECAY(J1)/DUMMY2
            !COORDS(J2,NP)=  COORDS(J2,NP)+  LOCALSTEP*RANDOMZ*EXPDECAY(J1)/DUMMY2
!!        ELSE IF (COOP.AND.(NATOMS-NCORE(NP).GE.MAX(2,NCOOP))) THEN
         !ELSE IF (COOP.AND.(NCORE(NP).GE.MAX(2,NCOOP))) THEN
!        ELSE IF (COOP) THEN
          !  J2=3*INDEXD(J1)
            !DUMMY=CDIST(J1)
            !IF (NCORE(NP).GT.0) THEN
               !J2=3*INDEXD(J1)
               !DUMMY=CDIST(J1)
            !ENDIF
!!           WRITE(LFH,'(A,2I6,F15.5,I6)') 'takestep> J1,J2/3,dist,NCOOPDONE=',J1,J2/3,DUMMY,NCOOPDONE
            !IF ((NCOOPDONE.LE.NCOOP+1).AND.(DUMMY.LT.COOPCUT)) THEN
               !NCOOPDONE=NCOOPDONE+1
               !WRITE(LFH,'(A,I6,A,F12.4)') 'takestep> coop move for atom ',J2/3,' cdist=',DUMMY
               !COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOMX
               !COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOMY
               !COORDS(J2,NP)=  COORDS(J2,NP)+  LOCALSTEP*RANDOMZ
!!
!! Then a random move as well!
!!
!!              COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
!!              COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
!!              COORDS(J2,NP)=  COORDS(J2,NP)+  LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
            !ELSE
               !COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
               !COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
               !COORDS(J2,NP)=  COORDS(J2,NP)+  LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
            !ENDIF
         ELSE IF (.NOT.FIXD) THEN
!
!  Tried changing to scale steps according to distance from CM. Maximum
!  allowed shift is linear with this distance. Worse.
!  Now try moving atoms according to how strongly bound they are.
!
           RANDOM=DPRAND()
           IF ((VMIN-VMAX.EQ.0.0D0).OR.(EFAC.EQ.0.0D0)) THEN
                 IF (NCORE(NP).GT.0) THEN ! project out radial component and rescale to same length
                                          ! Only for shell atoms, not core, if they are moved at all.
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
                    COORDS(J2-2,NP)=COORDS(J2-2,NP)+DELX
                    COORDS(J2-1,NP)=COORDS(J2-1,NP)+DELY
                    COORDS(J2,NP)=COORDS(J2,NP)+DELZ
                 ELSE
!
!	csw34> MAKE BASIC RANDOM CARTESIAN DISPLACEMENT MOVE
!              I've added the below IF line so that when MOVABLEATOMS is specified, cartesian displacements
!              specified by STEP will only be applied to atoms in that the 'movableatoms' file. It could be
!              that this needs changing into a new keyword to allow you to rotate bits and apply cartesian
!              moves to the whole system - something to keep in mind
!
                 RANDOM=(DPRAND()-0.5D0)*2.0D0
                 IF (GAUSST) RANDOM=RANDOM/GKSMALL(J2-2) ! scale gauss steps by 1/GKSMALL
                 COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM
                 RANDOM=(DPRAND()-0.5D0)*2.0D0
                 IF (GAUSST) RANDOM=RANDOM/GKSMALL(J2-1)
                 COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM
                 RANDOM=(DPRAND()-0.5D0)*2.0D0
                 IF (GAUSST) RANDOM=RANDOM/GKSMALL(J2)
                 COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM
!
!	csw34> END OF BASIC RANDOM CARTESIAN DISPLACEMENT MOVE
!
              ENDIF
           ELSE 
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              DUMMY=1.0D0+EAMP*TANH(-2.0D0*EFAC*(VAT(J1,NP)-(VMAX+VMIN)/2.0D0)/(VMAX-VMIN))
              COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM*DUMMY
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM*DUMMY
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM*DUMMY
           ENDIF
         ENDIF
!13       CONTINUE
      ENDDO

!
!  Preserve centre of mass if required.
!
      IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,NP))
      IF (FIXCOM.AND.(.NOT.SEEDT)) CALL CENTRECOM(COORDS(:3*NATOMS,NP))

      RETURN
      END


