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
MODULE NEWNEBMODULE
     IMPLICIT NONE
     CONTAINS

     SUBROUTINE NEWNEB(REDOPATH,TSREDO,EINITIAL,QQ,EFINAL,FINFIN,MOREP,NATOMSIN,NOPTIN,NINTSIN)
          USE PORFUNCS
          USE NEBDATA
          USE KEYNEB
          USE MINIMISER1
          USE MINIMISER2
          USE MINIMISER3
          USE NEBOUTPUT
          USE NEBUTILS
          USE KEY, ONLY : UNRST, GROWSTRINGT, FREEZENODEST, DESMDEBUG, &
               & NEBMUPDATE, MUPDATE, DEBUG, BFGSSTEPS, NEBRESEEDT, &
               & INTCONMAX, ORDERI, ORDERJ, EPSALPHA, REDOBFGSSTEPS, & 
               & NREPMAX, DISTREF, NEBKINITIAL, ADDREPT, REPPOW, REDOTSIM, MIN1REDO, MIN2REDO, PUSHOFF, &
               & CONI, CONJ, CONDISTREF, AMHT, NUMGLY, REPI, REPJ, BULKT, D1INIT, D2INIT, &
               & REDOKADD, REDOPATH1, &
               & REDOPATH2, NREPI, NREPJ, REPCUT, NREPCUT
          USE GROWSTRINGUTILS, ONLY: GROWSTRING, TOTSTEPS
          USE GSDATA, ONLY : KEYGSPRINT
          USE MODGUESS,ONLY: GUESSPATHT,NINTERP
          USE MODMEC,ONLY: MECCANOT          
          USE INTCOMMONS, ONLY : DESMINT, INTINTERPT, NINTIM, NDIH, DIHINFO, ALIGNDIR, PREVDIH, NINTC
          USE INTCUTILS, ONLY : INTINTERPOLATE, CART2INT
          USE SPFUNCTS, ONLY : DUMPCOORDS
          USE NEBTOCONNECT
          USE AMHGLOBALS, ONLY : NMRES
          USE COMMONS,ONLY: PARAM1,PARAM2,PARAM3,REDOPATHNEB,ZSYM
          IMPLICIT NONE

          COMMON /OLDC/ EMAX

          DOUBLE PRECISION,INTENT(IN)           :: EINITIAL, EFINAL
          DOUBLE PRECISION,DIMENSION(:)         :: QQ,FINFIN
          LOGICAL,INTENT(IN),OPTIONAL  :: MOREP
          INTEGER,INTENT(IN),OPTIONAL  :: NATOMSIN,NOPTIN,NINTSIN

          INTEGER :: J1,JMAX, NPERSIST, ITDONE, K, I, J2
          DOUBLE PRECISION :: EMAX
          DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: MYPTS ! JMC
          DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: VNEW, LCOORDS ! JMC
          LOGICAL PERSISTENT(NIMAGE+2)
          LOGICAL REDOPATH, MFLAG, PTEST, LPTEST, LRESET
          DOUBLE PRECISION ENERGY, RMS2, EREAL, TSREDO(*)

          ! efk: for growstrings and internals
          DOUBLE PRECISION, ALLOCATABLE :: DELTAX(:)
          DOUBLE PRECISION, POINTER :: TANPTR(:,:)
          LOGICAL :: GSMFLAG
          LOGICAL :: FAILED

          LOGICAL :: KNOWE, KNOWG, KNOWH 
          COMMON /KNOWN/ KNOWE, KNOWG, KNOWH 

          IF (DESMDEBUG) THEN
          ! output coordinates of endpoints we're trying to connect          
             CALL DUMPCOORDS(QQ,'tryconnect.A.xyz', .FALSE.)
             CALL DUMPCOORDS(FINFIN,'tryconnect.B.xyz', .FALSE.)
          ENDIF
          
          CALL MYCPU_TIME(STARTTIME,.TRUE.)
          ! setup parameters
          ! Natoms,Nopt,Nints,Nimage
          IF (PRESENT(NATOMSIN)) THEN
               NATOMS=NATOMSIN
          ELSE
               NATOMS=SIZE(QQ)/3
          ENDIF
          IF (NATOMS<=0) THEN
               PRINT '(1x,a)', 'Number of atoms is less or equal to zero. Stop.'
               CALL TSUMMARY
               STOP
          ELSE IF (DEBUG) THEN
               PRINT *, 'newneb> Number of atoms = ',Natoms
          ENDIF
          ALLOCATE(ORDERI(NREPMAX),ORDERJ(NREPMAX),EPSALPHA(NREPMAX),DISTREF(NREPMAX),REPPOW(NREPMAX))
          ALLOCATE(BADIMAGE(NIMAGE+2),BADPEPTIDE(NIMAGE+2))
          ADDREPT=.FALSE.
          IF (PRESENT(NOPTIN)) THEN
               NOPT=NOPTIN
          ELSE
             IF (DESMINT) THEN
                NOPT = NINTC
             ELSE IF (AMHT) THEN
                NOPT = 3*(NMRES*3)-NUMGLY*3 
             ELSE
                NOPT=3*NATOMS
             ENDIF
          ENDIF
          IF (PRESENT(NINTSIN)) THEN
               NINTS=NINTSIN
          ENDIF
          IF (NIMAGE<=0) THEN
               PRINT '(1x,a)', 'Number of images is less or equal to zero. Stop.'
               CALL TSUMMARY
               STOP
          ENDIF
          ! printing
          IF (PRESENT(MOREP)) THEN
               MOREPRINTING=MOREP
          ENDIF
          IF (DEBUG.OR.DESMDEBUG) MOREPRINTING=.TRUE.
          IF (MOREPRINTING) THEN
             IF (GROWSTRINGT) THEN
                CALL KEYGSPRINT(.FALSE.)
             ELSE
                CALL ALLKEYNEBPRINT
             ENDIF
             PRINT*
          ENDIF

          IF (UNRST) THEN
               ALLOCATE(MYPTS(3*NATOMS*NIMAGE)) ! JMC
               GRADTYPE="dnebu"
               TANTYPE=4
          ENDIF
          IF (OLDCONNECT) OPTIMIZETS = .FALSE.
          BADTAU=.FALSE.
          
          ! set up arrays
          ALLOCATE(XYZ(NOPT*(NIMAGE+2)),GGG(NOPT*(NIMAGE+2)),SSS(NOPT*(NIMAGE+2)),EEE(NIMAGE+2), &
   &               RRR(NIMAGE+2),TANVEC(NOPT,NIMAGE),DVEC(NIMAGE+1),NEWNEBK(NIMAGE+1),DEVIATION(NIMAGE+1),STEPIMAGE(NIMAGE))

!         NEWNEBK(1:NIMAGE+1)=NEBK
          NEWNEBK(1:NIMAGE+1)=NEBKINITIAL

          IF (DESMINT) THEN
             ALLOCATE(XYZCART(3*NATOMS*(NIMAGE+2)), GGGCART(3*NATOMS*(NIMAGE+2)), TRUEGRAD(3*NATOMS*(NIMAGE+2)))
             ALLOCATE(DIHINFO(NIMAGE+2,NDIH))
             XCART => XYZCART(3*NATOMS+1:3*NATOMS*(NIMAGE+1))
             GCART => GGGCART(3*NATOMS+1:3*NATOMS*(NIMAGE+1))
             DIHINFO(:,:) = 0.0D0
          ELSE
             ALLOCATE(TRUEGRAD(NOPT*(NIMAGE+2)))
          ENDIF          

          X         => XYZ(NOPT+1:NOPT*(NIMAGE+1))
          EIMAGE    => EEE(2:NIMAGE+1)
          G         => GGG(NOPT+1:NOPT*(NIMAGE+1))
          GSPR      => SSS(NOPT+1:NOPT*(NIMAGE+1))
          RMSFIMAGE => RRR(2:NIMAGE+1)
          TANPTR => TANVEC
          EEE = 0.0D0
          EEE(1)=EINITIAL
          EEE(NIMAGE+2)=EFINAL
          IF (DESMINT.AND..NOT.GROWSTRINGT) THEN
          PRINT*, "newneb>"
             XYZCART(:3*NATOMS) = QQ
             XYZCART(3*NATOMS*(NIMAGE+1)+1:) = FINFIN

             PREVDIH => DIHINFO(1,:)
             CALL CART2INT(QQ,XYZ(:NOPT))

             DO J1 = 2,NIMAGE+2
                ! align all other dihedrals to start
                DIHINFO(J1,:) = DIHINFO(1,:)
             ENDDO

             ALIGNDIR = .TRUE.
             PREVDIH => DIHINFO(NIMAGE+2,:)
             CALL CART2INT(FINFIN,XYZ(NOPT*(NIMAGE+1)+1:))
             ALIGNDIR = .FALSE.
          ELSE
             XYZ(:NOPT)=QQ
             XYZ(NOPT*(NIMAGE+1)+1:)=FINFIN
          ENDIF
          TANPTR => TANVEC

          IF (FREEZENODEST.OR.NEBRESEEDT) THEN
             ALLOCATE(IMGFREEZE(NIMAGE))
             IMGFREEZE(:) = .FALSE.
          ENDIF

          IF(GROWSTRINGT) THEN
             IF((DESMINT.AND.NOPT.NE.NINTC).OR.(.NOT.DESMINT.AND.NOPT.NE.3*NATOMS)) THEN
                print*, 'NOPT must be equal to 3*NATOMS or NINTC to use growstring.'
                print*, DESMINT, NINTC, 3*NATOMS, NOPT
                STOP
             ENDIF
             IF (DESMINT) THEN
                CALL GROWSTRING(QQ, FINFIN, NIMAGE, XCART, EIMAGE, TANPTR,RMS,GSMFLAG)
                XYZ(:) = 0.0D0
             ELSE
                CALL GROWSTRING(QQ, FINFIN, NIMAGE, X, EIMAGE, TANPTR,RMS,GSMFLAG)
                CALL IMAGEDISTRIBUTION
             ENDIF
             NITERDONE = TOTSTEPS
             
          ELSE
             ! construct the band
             IF (READGUESS.OR.(GUESSPATHT.AND.UNRST.AND.(NINTERP.GT.1)).OR.(MECCANOT)) THEN
                CALL RWG("r",.True.,1)
                READGUESS = .FALSE.
             ELSE
                IF (UNRST) THEN ! JMC
                   CALL UNRESDIHENEB(QQ,FINFIN,MYPTS)
                   XYZ(NOPT+1:NOPT*(NIMAGE+1))=MYPTS(1:NOPT*NIMAGE)
                ELSEIF (REDOPATHNEB) THEN
                   REDOKADD=.TRUE.
                   REDOPATH1=.TRUE.
                   REDOTSIM=NIMAGE*D1INIT/(D1INIT+D2INIT)+1
                   XYZ(NOPT*REDOTSIM+1:NOPT*(REDOTSIM+1))=TSREDO(1:NOPT)
                   ALLOCATE(DELTAX(NOPT))
                   IF (BULKT) THEN
                      DO K=1,NATOMS
                         DELTAX(3*(K-1)+1)=XYZ(NOPT*REDOTSIM+3*(K-1)+1) - XYZ(3*(K-1)+1) &
  &                          -PARAM1*NINT((XYZ(NOPT*REDOTSIM+3*(K-1)+1) - XYZ(3*(K-1)+1))/PARAM1)
                         DELTAX(3*(K-1)+2)=XYZ(NOPT*REDOTSIM+3*(K-1)+2) - XYZ(3*(K-1)+2) &
  &                          -PARAM2*NINT((XYZ(NOPT*REDOTSIM+3*(K-1)+2) - XYZ(3*(K-1)+2))/PARAM2)
                         DELTAX(3*(K-1)+3)=XYZ(NOPT*REDOTSIM+3*(K-1)+3) - XYZ(3*(K-1)+3) &
  &                          -PARAM3*NINT((XYZ(NOPT*REDOTSIM+3*(K-1)+3) - XYZ(3*(K-1)+3))/PARAM3)
                      ENDDO
                      DELTAX(1:NOPT)=DELTAX(1:NOPT)/REDOTSIM
                   ELSE
                      DELTAX(1:NOPT) = ( XYZ(NOPT*REDOTSIM+1:NOPT*(REDOTSIM+1) ) - XYZ(1:NOPT) )/REDOTSIM
                   ENDIF
                   DO I=2,REDOTSIM
                      XYZ(NOPT*(I-1)+1:NOPT*I) = XYZ(1:NOPT) + DELTAX*(I-1)
                   ENDDO
                   IF (.NOT.ALLOCATED(VNEW)) ALLOCATE(VNEW(NOPT))
                   IF (.NOT.ALLOCATED(LCOORDS)) ALLOCATE(LCOORDS(NOPT))
                   IF (DEBUG) PRINT '(A)',' newneb> minimising on the start side to make DNEB images'
                   LPTEST=.TRUE.
                   LCOORDS(1:NOPT)=TSREDO(1:NOPT)+PUSHOFF*(MIN1REDO(1:NOPT)-TSREDO(1:NOPT))/D1INIT
                   KNOWE=.FALSE.; KNOWG=.FALSE.
                   LRESET=.TRUE.
                   DO I=REDOTSIM,2,-1
                      CALL MYLBFGS(NOPT,MUPDATE,LCOORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,REDOBFGSSTEPS,LRESET, &
  &                                ITDONE,LPTEST,VNEW,.TRUE.,.FALSE.)
                      LRESET=.FALSE.
                      XYZ(NOPT*(I-1)+1:NOPT*I)=LCOORDS(1:NOPT)
                      IF (DEBUG) PRINT '(A,I6,A,F20.10)',' newneb> energy of image ',I,' is ',EREAL
                   ENDDO
                   REDOPATH1=.FALSE.
                   REDOPATH2=.TRUE.
                   IF (BULKT) THEN
                      DO K=1,NATOMS
                         DELTAX(3*(K-1)+1)=XYZ(NOPT*(NIMAGE+1)+3*(K-1)+1) - XYZ(NOPT*REDOTSIM+3*(K-1)+1) &
  &                          -PARAM1*NINT((XYZ(NOPT*(NIMAGE+1)+3*(K-1)+1) - XYZ(NOPT*REDOTSIM+3*(K-1)+1))/PARAM1)
                         DELTAX(3*(K-1)+2)=XYZ(NOPT*(NIMAGE+1)+3*(K-1)+2) - XYZ(NOPT*REDOTSIM+3*(K-1)+2) &
  &                          -PARAM2*NINT((XYZ(NOPT*(NIMAGE+1)+3*(K-1)+2) - XYZ(NOPT*REDOTSIM+3*(K-1)+2))/PARAM2)
                         DELTAX(3*(K-1)+3)=XYZ(NOPT*(NIMAGE+1)+3*(K-1)+3) - XYZ(NOPT*REDOTSIM+3*(K-1)+3) &
  &                          -PARAM3*NINT((XYZ(NOPT*(NIMAGE+1)+3*(K-1)+3) - XYZ(NOPT*REDOTSIM+3*(K-1)+3))/PARAM3)
                      ENDDO
                      DELTAX(1:NOPT)=DELTAX(1:NOPT)/(NIMAGE-REDOTSIM+1)
                   ELSE
                      DELTAX(1:NOPT) = ( XYZ(NOPT*(NIMAGE+1)+1:NOPT*(NIMAGE+2))    &
  &                                    - XYZ(NOPT*REDOTSIM+1:NOPT*(REDOTSIM+1)) )/(NIMAGE-REDOTSIM+1)
                   ENDIF
                   DO I=REDOTSIM+2,NIMAGE+1
                      XYZ(NOPT*(I-1)+1:NOPT*I) = TSREDO(1:NOPT) + DELTAX*(I-REDOTSIM-1)
                   ENDDO
                   IF (DEBUG) PRINT '(A)',' newneb> minimising on the finish side to make DNEB images'
                   LCOORDS(1:NOPT)=TSREDO(1:NOPT)+PUSHOFF*(MIN2REDO(1:NOPT)-TSREDO(1:NOPT))/D2INIT
                   KNOWE=.FALSE.; KNOWG=.FALSE.
                   LRESET=.TRUE.
                   DO I=REDOTSIM+2,NIMAGE+1
                      CALL MYLBFGS(NOPT,MUPDATE,LCOORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,REDOBFGSSTEPS,LRESET, &
  &                                ITDONE,LPTEST,VNEW,.TRUE.,.FALSE.)
                      LRESET=.FALSE.
                      XYZ(NOPT*(I-1)+1:NOPT*I)=LCOORDS(1:NOPT)
                      IF (DEBUG) PRINT '(A,I6,A,F20.10)',' newneb> energy of image ',I,' is ',EREAL
                   ENDDO
                   DEALLOCATE(DELTAX,VNEW,LCOORDS)
                   KNOWE=.FALSE.; KNOWG=.FALSE.; KNOWH=.FALSE.
                   REDOKADD=.FALSE.
                   REDOPATH2=.FALSE.
                ELSEIF (REDOPATH) THEN
                   XYZ(NOPT+1:NOPT*2) = TSREDO(1:NOPT)
                   EEE(2)=1.0D100
                ELSE
                   IF (INTINTERPT) THEN
                      CALL INTINTERPOLATE(QQ,FINFIN,NINTIM,NIMAGE,X,DESMDEBUG,FAILED)
                   ELSE
                      CALL MAKEIMAGE(EINITIAL,EFINAL,QQ,FINFIN)
!                      PRINT *, 'XYZ'
!                      DO J1 = 1, NIMAGE+2
!                         DO J2 = 3,3 !1, NATOMS/2
!                            K = (J1-1)*3*NATOMS+J2*3
!                            PRINT *, K, XYZ(K)
!                         ENDDO
!                      ENDDO
!                      STOP
                   ENDIF
                ENDIF
             ENDIF
             IF (UNRST) DEALLOCATE(MYPTS) ! JMC

             ! preoptimise if requested        
             IF (SQVVGUESS) THEN
                STEPTOT = 5.0D0 ! TO AVOID GRADIENT SCALING DURING SQVV
                CALL NEBSQVV(NOPT*NIMAGE)
             ENDIF
             
             NPERSIST=0
             IF ((.NOT.REDOPATH).OR.REDOPATHNEB) THEN
                SELECT CASE(MINTYPE)
                CASE("lbfgs")
                   IF (UNRST) THEN
                      CALL NEBBFGSINT(NINTS*NIMAGE,NEBMUPDATE)
                   ELSE
                      CALL NEBBFGS(NOPT*NIMAGE,NEBMUPDATE,NPERSIST,PERSISTENT)
                   END IF
                CASE("sqvv")
                   CALL NEBSQVV(NOPT*NIMAGE)
                END SELECT
             ENDIF
          ENDIF
          
          ! save final NEB coordinates and energy profile
          IF (DEBUG) THEN
             PRINT '(A,F12.4)',' newneb> mean image separation is ',SEPARATION/(NIMAGE+1)
!            DO J1=1,NIMAGE+1
!               PRINT '(A,F12.4,A,I8,A,F12.4)',' newneb> NEB k is ',NEWNEBK(J1),' for gap ',J1,' value=',DVEC(J1)
!            ENDDO
          ENDIF

          IF (DUMPNEBEOS) CALL WRITEPROFILE(0)
!          PRINT *, 'XYZ'
!          DO J1 = 1, NIMAGE+2
!             DO J2 = 3,3 !1, NATOMS/2
!                K = (J1-1)*3*NATOMS+J2*3
!                PRINT *, K, XYZ(K)
!             ENDDO
!          ENDDO
!          STOP
          IF (DUMPNEBXYZ) CALL RWG("w",.False.,0)
          IF (DUMPNEBPTS) CALL SAVEBANDCOORD

          IF (OLDCONNECT) THEN ! FIND THE HIGHEST ENERGY IMAGE
             IF (DESMINT) THEN
                print*, 'newneb>> ERROR! OLDCONNECT not implemented with DESMINT'
                STOP
             ENDIF
               EMAX=MAXVAL(EIMAGE)
               DO J1=1,NIMAGE
                    IF (EMAX == EIMAGE(J1)) JMAX=J1
               ENDDO
               QQ = X(NOPT*(JMAX-1)+1:NOPT*JMAX)
          ENDIF
          
          CALL MYCPU_TIME(ENDTIME,.FALSE.)

          NMINFOUND=0
          NTSFOUND=0
!
!  Current structure precludes searching the NEB profile for
!  both minima and ts. Could perhaps change this. Current philosophy is
!  that if we have persistent minima we should start the whole DNEB again
!  without ts searches.
!
          IF (NPERSIST.GT.0) THEN
             PTEST=.FALSE.
             IF (DEBUG) PTEST=.TRUE.
             PRINT '(A,I8)',' newneb> number of persistent minima in DNEB profile=',NPERSIST 
             DO J1=2,NIMAGE+1
                IF (PERSISTENT(J1)) THEN
                   KNOWG=.FALSE.
                   KNOWE=.FALSE. ! could use EEE value
                   IF (.NOT.ALLOCATED(VNEW)) ALLOCATE(VNEW(NOPT))
                   PRINT '(A,I8,A,F20.10)',' newneb> minimising image ',J1,' initial energy=',EEE(J1)
                   CALL MYLBFGS(NOPT,MUPDATE,XYZ(NOPT*(J1-1)+1:NOPT*J1),.FALSE., &
   &                            MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,.TRUE.,ITDONE,PTEST,VNEW,.TRUE.,.FALSE.)
                   IF (MFLAG) THEN
                      NMINFOUND=NMINFOUND+1
!
!  We have to communicate the minima found back to tryconnect using the data structure
!  set up for new transition states. 
!  Added new variable MINFOUND to allow for this check in tryconnect.
!  It seems impossible to make newneb see isnewmin and addnewmin for some reason.
!
                      ALLOCATE(MINFOUND(NMINFOUND)%E,MINFOUND(NMINFOUND)%COORD(NOPT))
                      MINFOUND(NMINFOUND)%COORD(1:NOPT)=XYZ(NOPT*(J1-1)+1:NOPT*J1)
                      MINFOUND(NMINFOUND)%E=EREAL
                      WRITE(987,'(I6)') NATOMS
                      WRITE(987,'(A,I5)') 'image ',J1
                      WRITE(987,'(A,3G20.10)') (ZSYM(J2), MINFOUND(NMINFOUND)%COORD(3*(J2-1)+1:3*(J2-1)+3),J2=1,NATOMS)
                   ENDIF
                   DEALLOCATE(VNEW)
                ENDIF
             ENDDO
          ELSE
             CALL PRINTSUMMARY
             CALL TSLOCATOR
          ENDIF

          NULLIFY(X,EIMAGE)
          DEALLOCATE(DVEC,NEWNEBK,XYZ,EEE,GGG,TRUEGRAD,SSS,RRR,DEVIATION,TANVEC,STEPIMAGE,ORDERI,ORDERJ,EPSALPHA,DISTREF,REPPOW)
          DEALLOCATE(BADIMAGE,BADPEPTIDE)
          IF (FREEZENODEST) DEALLOCATE(IMGFREEZE)
          IF (DESMINT) THEN
             NULLIFY(XCART, GCART)
             DEALLOCATE(XYZCART,GGGCART,DIHINFO)
          ENDIF

     END SUBROUTINE NEWNEB
END MODULE NEWNEBMODULE
