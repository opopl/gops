!  GMIN: A program for finding global minima
!  Copyright (C) 1999-2006 David J. Wales
!  This file is part of GMIN.
!
!  GMIN is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  GMIN is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
! 
!---======================================---
      SUBROUTINE PTBASINSAMPLING

      USE MODCHARMM
      USE COMMONS
      USE tetherfunc

      IMPLICIT NONE

#ifdef MPI
!op226> MPI section {{{ 
      INCLUDE 'mpif.h'
!op226 Variable declarations      {{{
      INTEGER, PARAMETER :: NHIST=100, NHISTE=1000
      INTEGER :: IACCEPT(0:NPAR-1), MPIERR, J,K, IS(MPI_STATUS_SIZE)
      INTEGER NHISTQ4(NHIST,0:NPAR-1), NDUMMY, NTYPEA, NTOT, NH, IMESG, IQE, IQ4, IQ6, J1, J2, J3, J4,
     1        NHISTQ6(NHIST,0:NPAR-1), NHISTQE(NHISTE, 0:NPAR-1), IENR, 
     2        NOUT(0:NPAR-1), ITRAJ, ITRAJO,NEACCEPT, RNDSEED, NUPDATE,
     3        CONVERGED,LBFGS_ITERATIONS, JD, BININDEX, MINIMANUMBER(HBINS,0:NPAR-1),
     4        NHISTALLQ(NHIST, NHIST, 0:NPAR-1), IBININDEX, IBININDEX2, LOWESTDIRECT(HBINS),
     5        NHISTGRAND(NHISTE, NHIST, NHIST, 0:NPAR-1),LBFGS_ITERATIONSO,LVISITS(NENRPER),PEVISITS2(NENRPER, HBINS, 0:NPAR-1),
     7        NOUTQBIN, NOUTPEBIN, LBFGS_ITERATIONSSAVE, HBINMAX, HBINMIN
      DOUBLE PRECISION I
      DOUBLE PRECISION QVISITS(HBINS, 0:NPAR-1), PEVISITS(NENRPER,0:NPAR-1), NACCEPTPT(0:NPAR-1)

      DOUBLE PRECISION V(NATOMS), VO(NATOMS), TEMPTRAJ(0:NPAR-1), H(0:NPAR-1), BETA(0:NPAR-1), 
     1        EAV(0:NPAR-1), EAV2(0:NPAR-1), Q(3,NATOMS), Q4AV(0:NPAR-1), Q4AV2(0:NPAR-1), 
     2        Q6AV(0:NPAR-1), Q6AV2(0:NPAR-1), VENR(NENRPER), 
     3        HINIT(0:NPAR-1), EPSAB, EPSBB, SIGAB, SIGBB, X(NATOMS), Y(NATOMS), Z(NATOMS), 
     4        CTE, T, VOLD,VNEW, POTEL, GRAD(3*NATOMS), Q4, Q6, RANDOM, DPRAND, Q4MAX, Q6MAX , 
     5	      DQ4, DQ6, DHISTE, ENUL, XO(NATOMS), YO(NATOMS), ZO(NATOMS), DDX, DDY, DDZ, DE, 
     6        W, WCOMP, WAC, E, ER, DBETA, DELTA, CV, FQ4, FQ6, RMAX, DDXN, DDYN, DDZN, R2, 
     7        DUMMY,BINLABEL(HBINS), VNEWSAVE, DIHEORDERPARAM,SASAORDERPARAM, PEINT, HISTINT, 
     8        DIHEORDERPARAM_AV(0:NPAR-1), SASAORDERPARAM_AV(0:NPAR-1), MINCOORDS(3*NATOMS,NPAR), INTERIMCOORDS(3*NATOMS,1),
     9        INTERIMCOORDS2(3*NATOMS,1), DIST, VOLDSAVE, LESAVE, DOSSTATS(MAXIT,2), DOSSTATSO(MAXIT,2), 
     A        DCOORDS(3*NATOMS), DOSSTATSDUM(MAXIT,2), SR3

       DOUBLE PRECISION DIHEORDERPARAM_MIN, DIHEORDERPARAM_MAX, SASAORDERPARAM_MIN, SASAORDERPARAM_MAX, 
     1          DDIHE, DSASA, RGYR_MIN, RGYR_MAX, EINT_MIN, EINT_MAX, DRGYR, DEINT,
     2          RGYR_AV(0:NPAR-1), EINT_AV(0:NPAR-1), RGYR_AV2(0:NPAR-1), EINT_AV2(0:NPAR-1),
     3          ORDERPARAM1, ORDERPARAM2, DISTANCE(HBINS, 0:NPAR-1),
     4          DISTANCE_AV(HBINS,0:NPAR-1), DISTANCEOLD, NORM_PJ(0:NPAR-1), PJ(HBINS,0:NPAR-1),
     5          VMINOLD, VMINNEW, DELTA_EAV(0:NPAR-1), CX, CY, CZ,
     6        DIHEORDERPARAM_AV2(0:NPAR-1), SASAORDERPARAM_AV2(0:NPAR-1), RGYR, EINT, VMINNEWSAVE

      CHARACTER (LEN =256)  FILENAME, FILENAME2,FILENAME3,FILENAME4,FILENAME5,FILENAME6,
     1                      FILENAME7,FILENAME8,FILENAME9, FILENAME10, FILENAME11, FILENAME12, 
     2                      FILENAME100,FILENAME101,FILENAME102, FILENAME103, FILENAME104, FILENAME105 
      CHARACTER (LEN= 80)  ISTR, SDUMMY
      CHARACTER (LEN= 80)  CFNAME
      LOGICAL EXCHANGE, EXCHANGEACCEPT, FITS, NEWENERGY, EVAP, EVAPREJECT, RECOUNT, YESNO

      DOUBLE PRECISION PREVSTEPS

      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      COMMON /MYPOT/ POTEL
      COMMON /EV/ EVAP, EVAPREJECT
! }}}

      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NDUMMY,MPIERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNODE,MPIERR)
      NOUTQBIN=0
      NOUTPEBIN=0

      WRITE(LFH, '(A,I10,A,I10)') "bspt> This is processor", mynode+1, " of", NPAR
      WRITE(LFH, '(A,I10)') 'bspt> Number of atoms', natoms
      IF (FIXSTEP(1)) FIXSTEP(2:NPAR)=.TRUE.
      IF (FIXTEMP(1)) FIXTEMP(2:NPAR)=.TRUE.
      IF (PERIODIC) THEN
         WRITE(LFH, '(A,6G20.10)') 'bspt> Binary data', ntypea, epsab, epsbb, sigab, sigbb, cutoff
         WRITE(LFH, '(A,3G20.10)') 'bspt> Box data', boxlx, boxly, boxlz
      ELSEIF(CHRMMT) THEN
         WRITE(LFH, '(A)') 'bspt> CHARMM job'
         IF (FIXSTEP(MYNODE+1).OR.FIXBOTH(MYNODE+1)) THEN
            WRITE(LFH, '(A,G20.10,A)') 'bspt> Starting stepsize ', STEP(MYNODE+1),' (fixed)'
         ELSE
            WRITE(LFH, '(A,G20.10,A)') 'bspt> Starting stepsize ', STEP(MYNODE+1),
     &                            ' (will be dynamically adjusted during equilibration)'
         ENDIF
      ELSEIF (.NOT.MODEL1T) THEN
         WRITE(LFH, '(A,2G20.10)') 'bspt> Radius and Radius**2: ',SQRT(RADIUS),radius
      ENDIF

      ITRAJ=MYNODE
      NEACCEPT=0
!     MINIMANUMBER=0 ! Number of distinct minima for quench bins in each replica. Not currently used.
      IENR=0 ! DJW APPEARED TO BE UNINITIALISED ?
      PEINT=(PTEMAX-PTEMIN)/NENRPER
      MAXEFALL=-1.0D100
      HISTINT=(HISTMAX-HISTMIN)/HBINS
      IF (MINDENSITYT.AND.(PTSTEPS*1.0D0.GT.0.0D0)) THEN
         WRITE(LFH, '(A,2I8)') 'bspt> WARNING **** PTSTEPS reset to zero for density of minima run' 
         PTSTEPS=0
      ENDIF
!
! Limit for quench energy above which steps will be rejected.
!
      HBINMAX=(MIN(BSPTQMAX,HISTMAX)-HISTMIN)/HISTINT+1
      HBINMIN=(MAX(BSPTQMIN,HISTMIN)-HISTMIN)/HISTINT
      WRITE(LFH, '(A,2I8)') 'bspt> Minimum and maximum quench bin limits: ',HBINMIN,HBINMAX
      DO J1=1, HBINS
!        BINLABEL(J1)=HISTMIN + HISTINT*(J1-0.5D0) ! these energies point to the middle of the quench bin
         BINLABEL(J1)=HISTMIN + HISTINT*(J1-1.0D0) ! these energies point to the bottom of the quench bin
      ENDDO

! Initialisation

      DO J1=1,NATOMS
         X(J1)=COORDS(3*(J1-1)+1,MYNODE+1)
         Y(J1)=COORDS(3*(J1-1)+2,MYNODE+1)
         Z(J1)=COORDS(3*(J1-1)+3,MYNODE+1)
      ENDDO
      DO J1=1,NATOMS
         Q(1,J1)=X(J1)
         Q(2,J1)=Y(J1)
         Q(3,J1)=Z(J1)
      ENDDO

      CTE=(LOG(PTTMAX/PTTMIN))/(NPAR-1)
      CTE=EXP(CTE)

      DO J1=0, NPAR-1
         TEMPTRAJ(J1)=PTTMIN*CTE**J1
         T=TEMPTRAJ(J1)
         BETA(J1)=1.0D0/T
      ENDDO
      IF (MYNODE.EQ.0) THEN
         OPEN(UNIT=10, FILE='temperatures',STATUS='UNKNOWN')
         WRITE(10,'(G20.10)') TEMPTRAJ(0:NPAR-1)
         CLOSE(10)
      ENDIF
!
! For restart we need to get the current configuration, its pe, the pe of the minimum it quenched to,
! if applicable, the number of steps already done, the maximum step size, and the Visits and Visits2
! histograms. If we dump using BSPTDUMPFRQ then we can restore from the last such file. We can work
! out what the last dump was once we know how many steps have been done!
!
      IF (BSPTRESTART) THEN
         WRITE (ISTR, '(I10)') MYNODE
         FILENAME12="bsptrestart."//TRIM(ADJUSTL(ISTR))
         OPEN(UNIT=1986+MYNODE,FILE=FILENAME12, STATUS="old", form="formatted")
         READ(1986+MYNODE,*) PREVSTEPS,VOLD,VMINOLD,STEP(MYNODE+1),
     &                       NACCEPTPT(MYNODE),NEACCEPT,NTOT,NOUTQBIN,NOUTPEBIN,NQ(MYNODE+1)
         WRITE (LFH,'(A,F20.1)') 'bspt> Number of previous steps=     ',PREVSTEPS
         WRITE (LFH,'(A,2F20.10)') 'bspt> Quench and instantaneous PE=',VMINOLD,VOLD
         WRITE (LFH,'(A,F20.10)') 'bspt> Step size=',STEP(MYNODE+1)
         WRITE (LFH,'(A,F15.1,2I10)') 'bspt> Accepted MC steps, PT steps and total PT steps=',NACCEPTPT(MYNODE),NEACCEPT,NTOT
         WRITE (LFH,'(A,3I10)') 'bspt> Quenches and instantaneous energies outside range=',NOUTQBIN,NOUTPEBIN
         WRITE (LFH,'(A,I10)') 'bspt> Total quenches=',NQ(MYNODE+1)
         DO J1=1,NATOMS
            J2=3*(J1-1)
            READ(1986+MYNODE,*) COORDS(J2+1,MYNODE+1), COORDS(J2+2,MYNODE+1), COORDS(J2+3,MYNODE+1)
         ENDDO
         CLOSE(1986+MYNODE)
         IF (BSPTDUMPFRQ.GT.0) THEN
            DUMMY=INT(PREVSTEPS/(1.0D0*BSPTDUMPFRQ))*1.0D0
!
!  Miscompiles unless we split up the calculation of DUMMY !!!
!
!           WRITE(LFH,*) 'bspt> PREVSTEPS,BSPTDUMPFRQ,DUMMY=',PREVSTEPS,BSPTDUMPFRQ,DUMMY
!           WRITE(LFH,*) 'bspt> INT,arg=',INT(PREVSTEPS/(1.0D0*BSPTDUMPFRQ)),PREVSTEPS/(1.0D0*BSPTDUMPFRQ)
            DUMMY=DUMMY*BSPTDUMPFRQ
!           WRITE(LFH,*) 'bspt> DUMMY=',DUMMY
            CALL FLUSH(LFH)
            WRITE (SDUMMY, '(F15.1)') DUMMY
            ISTR=TRIM(ADJUSTL(ISTR)) // '.' // TRIM(ADJUSTL(SDUMMY))
         ELSE
            ISTR=TRIM(ADJUSTL(ISTR)) ! this should be the final and only Visits file
         ENDIF
         IF (PREVSTEPS.GT.NEQUIL) THEN
            FILENAME101="Visits.his."//TRIM(ADJUSTL(ISTR))
            WRITE(LFH,'(A,A)') 'bspt> Reading restart information from ',TRIM(ADJUSTL(FILENAME101))
            CALL FLUSH(LFH)
            OPEN(UNIT=1986+MYNODE,FILE=FILENAME101, STATUS="unknown", form="formatted")
            READ(1986+MYNODE, '(G20.10)') DUMMY
            READ(1986+MYNODE, '(A)') SDUMMY
            DO K=1, NENRPER
               READ(1986+MYNODE,*) DUMMY,PEVISITS(K,MYNODE)
            ENDDO
            READ(1986+MYNODE, '(A)') SDUMMY
            DO K=1, HBINS
               WRITE(1986+MYNODE,*) DUMMY, QVISITS(K,MYNODE)
            ENDDO
            CLOSE(1986+MYNODE)
         ENDIF
         IF (BSPT.AND.(PREVSTEPS.GT.NEQUIL+PTSTEPS).AND.(.NOT.MINDENSITYT)) THEN
            FILENAME101="Visits2.his."//TRIM(ADJUSTL(ISTR))
            WRITE(LFH,'(A,A)') 'bspt> Reading restart information from ',TRIM(ADJUSTL(FILENAME101))
            OPEN(UNIT=1986+MYNODE,FILE=FILENAME101,STATUS="unknown",FORM="UNFORMATTED")
            READ(1986+MYNODE) PEVISITS2(1:NENRPER,1:HBINS,MYNODE)
            CLOSE(1986+MYNODE) 
         ENDIF
         CONVERGED=1
      ELSE
12       CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
         VOLD=POTEL
         WRITE(LFH,'(2(A,G20.10))') 'bspt> Initial configuration energy is ',VOLD
!        WRITE(LFH,'(2(A,G20.10))') 'bspt> NOT QUENCHING initial config'
         CALL QUENCH(.FALSE.,MYNODE+1,LBFGS_ITERATIONS,DUMMY,NDUMMY,CONVERGED,COORDS(:,MYNODE+1),DOSSTATS)
         WRITE(LFH,'(A,G20.10)') 'bspt> Initial configuration has been quenched to energy ',POTEL
!
!  If we start with an invalid configuration we may never recover!
!
         IF (.NOT.(CHRMMT.OR.MODEL1T)) THEN
            DO K=1,NATOMS
               DIST=COORDS(3*(K-1)+1,MYNODE+1)**2+COORDS(3*(K-1)+2,MYNODE+1)**2+COORDS(3*(K-1)+3,MYNODE+1)**2
               IF (DIST.GT.RADIUS) THEN
                  IF (MOD(I-1.0D0,1.0D0*PRTFRQ).EQ.0) WRITE(LFH,'(A,I6,A)') 'bspt> Atom ',K,
     &                       ' outside container, reseed and try again'
                  SR3=DSQRT(3.0D0)
                  DO J1=1,NATOMS
                     RANDOM=(DPRAND()-0.5D0)*2.0D0
                     COORDS(3*(K-1)+1,MYNODE+1)=RANDOM*DSQRT(RADIUS)/SR3
                     RANDOM=(DPRAND()-0.5D0)*2.0D0
                     COORDS(3*(K-1)+2,MYNODE+1)=RANDOM*DSQRT(RADIUS)/SR3
                     RANDOM=(DPRAND()-0.5D0)*2.0D0
                     COORDS(3*(K-1)+3,MYNODE+1)=RANDOM*DSQRT(RADIUS)/SR3
                  ENDDO
                  GOTO 12
               ENDIF
            ENDDO 
         ENDIF
      ENDIF

      LBFGS_ITERATIONS=0
      LBFGS_ITERATIONSO=0 ! so that we don't use quench DoS statistics from arbitrary high energy
      VOLD=POTEL
      DO J1=1,NATOMS
         X(J1)=COORDS(3*(J1-1)+1,MYNODE+1)
         Y(J1)=COORDS(3*(J1-1)+2,MYNODE+1)
         Z(J1)=COORDS(3*(J1-1)+3,MYNODE+1)
      ENDDO
      DO J1=1,NATOMS
         Q(1,J1)=X(J1)
         Q(2,J1)=Y(J1)
         Q(3,J1)=Z(J1)
      ENDDO

      IF (PERIODIC) THEN
         CALL QORDER_BLJ(Q,Q4,Q6)
      ELSE IF (CHRMMT) THEN
           IF (ODIHET) CALL CHCALCDIHE(DIHEORDERPARAM,COORDS(1:3*NATOMS,MYNODE+1))
           IF (OSASAT) CALL ORDER_SASA(SASAORDERPARAM,RPRO,COORDS(1:3*NATOMS:3,MYNODE+1),
     1             COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1))
           IF (ORGYT) CALL CHCALCRGYR(RGYR,COORDS(1:3*NATOMS:3,MYNODE+1),
     1             COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1),.FALSE.) 
           IF (OEINTT) CALL CHCALCEINT(EINT,COORDS(1:3*NATOMS:3,MYNODE+1),
     1             COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1),POTEL)
      ELSE
         CALL QORDER_LJ(Q,Q4,Q6)
      ENDIF

      WRITE(LFH, '(A,2G20.10)') 'bspt> Temperature range', TEMPTRAJ(0), TEMPTRAJ(NPAR-1)
      WRITE(LFH, '(A,G20.10)') 'bspt> This temperature trajectory=', TEMPTRAJ(MYNODE)
      IF (.NOT.BSPTRESTART) VMINOLD=VOLD
      IF ((.NOT.BSPTRESTART).AND.(NEQUIL+PTSTEPS.GT.0)) VMINOLD=0.0D0 
      WRITE (LFH,'(A,F20.1)') 'bspt> Number of equilibration steps=',NEQUIL
      WRITE (LFH,'(A,F20.1)') 'bspt> Number of PT steps=           ',PTSTEPS
      WRITE (LFH,'(A,F20.1)') 'bspt> Number of BSPT steps=         ',NQUENCH
      IF (MINDENSITYT) WRITE (LFH,'(A,F20.1)') 'bspt> Accept/reject based on quench energy'
      IF (PERIODIC) WRITE(LFH, '(A, 2G20.10)') 'bspt> Starting Q4, Q6=', Q4, Q6
      IF (CHRMMT) THEN
        IF (OSASAT)  WRITE(LFH, '(A, 2G20.10)') 'bspt> Starting Met-enk order params=', 
     1                                         DIHEORDERPARAM, SASAORDERPARAM 
        IF (ORGYT)  WRITE(LFH, '(A, G20.10)') 'bspt> Starting Rgyr order params=', 
     1                                          RGYR 
        IF (OEINTT) WRITE(LFH, '(A, G20.10)') 'bspt> Starting Eint order params=',
     1                                          EINT
      ENDIF

! Initialisation complete

      RNDSEED=2002+MYNODE
      CALL SDPRND(RNDSEED)
      RANDOM=DPRAND()
      WRITE(LFH, '(A, G20.10)') 'bspt> Starting random number=', RANDOM
      IF (FIXSTEP(MYNODE+1).OR.FIXBOTH(MYNODE+1)) THEN
         WRITE(LFH, '(A,G20.10,A)') 'bspt> Starting stepsize ', STEP(MYNODE+1),' (fixed)'
      ELSE
         WRITE(LFH, '(A,G20.10,A)') 'bspt> Starting stepsize ', STEP(MYNODE+1),' (will be dynamically adjusted)'
      ENDIF
      CALL FLUSH(LFH)

      IF (.NOT.BSPTRESTART) NTOT=0
      IF (.NOT.BSPTRESTART) NACCEPTPT(MYNODE)=0.0D0
      IF (.NOT.BSPTRESTART) QVISITS(:,MYNODE)=0.0D0
      IF (.NOT.BSPTRESTART) PEVISITS(:,MYNODE)=0.0D0
      IF (.NOT.BSPTRESTART) PEVISITS2(:,:,MYNODE)=0
      IACCEPT(MYNODE)=0
      EAV(MYNODE)=0.
      EAV2(MYNODE)=0.
      Q4AV(MYNODE)=0.
      Q6AV(MYNODE)=0.
      Q4AV2(MYNODE)=0.
      Q6AV2(MYNODE)=0.
      DIHEORDERPARAM_AV(MYNODE)=0.
      SASAORDERPARAM_AV(MYNODE)=0.
      DIHEORDERPARAM_AV2(MYNODE)=0.
      SASAORDERPARAM_AV2(MYNODE)=0.
      RGYR_AV(MYNODE)=0.
      EINT_AV(MYNODE)=0.
      RGYR_AV2(MYNODE)=0.
      EINT_AV2(MYNODE)=0.
      DO J1=1,NHIST
         NHISTQ4(J1,MYNODE)=0
         NHISTQ6(J1,MYNODE)=0
      ENDDO
      DO J1=1,NHISTE
         NHISTQE(J1,MYNODE)=0
      ENDDO
      NHISTALLQ(:,:,MYNODE)=0
      NHISTGRAND(:,:,:,MYNODE)=0
      DISTANCE(:,MYNODE)=0.0D0
      DISTANCE_AV(:,MYNODE)=0.0D0
      DELTA_EAV(MYNODE)=0.0D0

      Q4MAX=0.1
      Q6MAX=0.5 ! WILL HAVE TO CHANGE AS IS RATHER SYSTEM-SPECIFIC
      DIHEORDERPARAM_MIN=0.4D0
      DIHEORDERPARAM_MAX=1.0D0
      SASAORDERPARAM_MIN=300.0D0
      SASAORDERPARAM_MAX=800.0D0
      RGYR_MIN=4.0D0
      RGYR_MAX=12.0D0
      EINT_MIN=-35.0
      EINT_MAX=10.0
      DDIHE=(DIHEORDERPARAM_MAX-DIHEORDERPARAM_MIN)/(NHIST-1)
      DSASA=(SASAORDERPARAM_MAX-SASAORDERPARAM_MIN)/(NHIST-1)
      DRGYR=(RGYR_MAX-RGYR_MIN)/(NHIST-1)
      DEINT=(EINT_MAX-EINT_MIN)/(NHIST-1)
      DQ4=Q4MAX/(NHIST-1)
      DQ6=Q6MAX/(NHIST-1)
      DHISTE=(PTEMAX-PTEMIN)/(NHISTE-1)


      ENUL=VOLD

      IF (MYNODE.EQ.0) THEN

         OPEN(UNIT=10, FILE='distributions.header', FORM='formatted')
         WRITE(10,*) NATOMS
         WRITE(10,*) NPAR, PTSTEPS/NENRPER, NENRPER
         WRITE(10,*) ENUL
         CLOSE(10)
      ENDIF 
!
! NUPDATE specifies the interval for dynamically altering the maximum step size.
! Only used if step size isn't fixed.
!
      NUPDATE=100
!
! Main loop over steps.
!
!     DO I=1.0D0,NEQUIL+PTSTEPS+NQUENCH

      I=1.0D0
      IF (BSPTRESTART) I=PREVSTEPS
      DO 
         I=I+1.0D0
         IF (I.GT.NEQUIL+PTSTEPS+NQUENCH) EXIT
         RECOUNT=.FALSE.
         DO K=1, NATOMS
            XO(K)=X(K)
            YO(K)=Y(K)
            ZO(K)=Z(K)
         ENDDO
!
! If a move is rejected, the result for the current configuration
! should be recounted. This allows for the unsymmetrical forward
! and backward move probabilities for an atom near the surface of
! the constraining sphere.
!             
         IF (CHRMMT) THEN
!           CALL TAKESTEPCH(MYNODE+1)
!           DO K=1,NATOMS
!              X(K)=COORDS(3*(K-1)+1,MYNODE+1)
!              Y(K)=COORDS(3*(K-1)+2,MYNODE+1)
!              Z(K)=COORDS(3*(K-1)+3,MYNODE+1)
!           ENDDO
            DO K=1,NATOMS
               RANDOM=DPRAND()
               X(K) = X(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)
               RANDOM=DPRAND()
               Y(K) = Y(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)
               RANDOM=DPRAND()
               Z(K) = Z(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)
            ENDDO
         ELSE
            DO K=1,NATOMS
               RANDOM=DPRAND()
               X(K) = X(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)  
               RANDOM=DPRAND()
               Y(K) = Y(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)  
               RANDOM=DPRAND()
               Z(K) = Z(K) + 2.0D0*RANDOM*STEP(MYNODE+1)-STEP(MYNODE+1)  
            ENDDO
         ENDIF
!
! Should probably worry about centre of coordinates versus centre of mass
! for heteroatomic systems.
!
         IF (CENT) THEN
            CX=0.0D0; CY=0.0D0; CZ=0.0D0
            DO K=1,NATOMS
               CX=CX+X(K); CY=CY+Y(K); CZ=CZ+Z(K)
            ENDDO
            CX=CX/NATOMS; CY=CY/NATOMS; CZ=CZ/NATOMS
            DO K=1,NATOMS
               X(K)=X(K)-CX; Y(K)=Y(K)-CY; Z(K)=Z(K)-CZ
            ENDDO
         ENDIF 
         IF (.NOT.(CHRMMT.OR.MODEL1T)) THEN
            cloop: DO K=1,NATOMS
               DIST=X(K)**2+Y(K)**2+Z(K)**2
               IF (DIST.GT.RADIUS) THEN
                  IF (MOD(I-1.0D0,PRTFRQ*1.0D0).EQ.0.0D0) WRITE(LFH,'(A,I6,A)') 'bspt> Perturbed atom ',K,
     &                       ' outside container, recount previous configuration results'
                  RECOUNT=.TRUE. ! The right way to deal with rejected steps!
                  EXIT cloop
               ENDIF
            ENDDO cloop
         ENDIF
                
         DO K=1,NATOMS
            COORDS(3*(K-1)+1,MYNODE+1)=X(K)
            COORDS(3*(K-1)+2,MYNODE+1)=Y(K)
            COORDS(3*(K-1)+3,MYNODE+1)=Z(K)
         ENDDO
!
! At this point all we have done is take a step. The perturbed coordinates are in both
! COORDS and X, Y, Z. The old coordinates are in XO, YO, ZO.
!
! New and old quench energies will be in        VMINNEW VMINOLD
! New and old instantaneous energies will be in    VNEW    VOLD
!
         CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
         VNEW=POTEL
         IF (.NOT.MINDENSITYT) THEN
            WCOMP=(VNEW-VOLD)*BETA(MYNODE) ! use difference in instantaneous energies
            W=MIN(1.0D0,EXP(-WCOMP))
            RANDOM=DPRAND()
            IF (RANDOM.GT.W) RECOUNT=.TRUE. ! RECOUNT is initialised to .FALSE. at the top of the loop
         ENDIF
!
! Quenching part if required.
! COORDSO saves the perturbed coordinates before the quench in order to calculate
! a quench distance. Should no longer be needed. COORDS are used as scratch for quenches.
!
         COORDSO(:,MYNODE+1)=COORDS(:,MYNODE+1) 
         IF ((RECOUNT.AND.BSPT).OR.(I.LE.NEQUIL+PTSTEPS)) THEN
            VMINNEW=0.0D0
            LBFGS_ITERATIONS=0
         ENDIF
         IF (MINDENSITYT.OR.(I.GT.NEQUIL+PTSTEPS)) THEN
            IF ((.NOT.RECOUNT).AND.BSPT) THEN 
               IF (MOD(I,1.0D0*QUENCHFRQ).EQ.1.0D0*0) THEN 
!                 WRITE(LFH, '(A)') 'bspt> calling quench'
                  IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
                  CALL QUENCH(.FALSE.,MYNODE+1,LBFGS_ITERATIONS,DUMMY,NDUMMY,CONVERGED,COORDS(:,MYNODE+1),DOSSTATS)
                  IF (CONVERGED.NE.1) WRITE(LFH, '(A)') 'bspt> WARNING - quench did not converge' 
                  VMINNEW=POTEL
                  NQ(MYNODE+1)=NQ(MYNODE+1)+1
                  MINCOORDS(:,MYNODE+1)=COORDS(:,MYNODE+1) ! MINCOORDS contains quench coords
               ELSE
                  VMINNEW=0.0D0
                  LBFGS_ITERATIONS=0
                  EVAPREJECT=.FALSE.
                  CONVERGED=1
                  POTEL=HUGE(1.0D0)
               ENDIF
            ELSEIF (RECOUNT.AND.BSPT.AND.(VMINOLD.EQ.0.0D0).AND.(MOD(I,1.0D0*QUENCHFRQ).EQ.0.0D0)) THEN
!
! We might not have quenched at the step to be recounted, so we have to do so here.
!
               IF (DEBUG) WRITE(LFH, '(A,G20.10,A)') 'bspt> recounting step for previous configuration with VMINOLD=',
     &                                     VMINOLD,' need to call quench'
               DO K=1,NATOMS
                 COORDS(3*(K-1)+1,MYNODE+1)=XO(K)
                 COORDS(3*(K-1)+2,MYNODE+1)=YO(K)
                 COORDS(3*(K-1)+3,MYNODE+1)=ZO(K)
               ENDDO
!              WRITE(LFH, '(A)') 'bspt> calling quench for recount'
               IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
               CALL QUENCH(.FALSE.,MYNODE+1,LBFGS_ITERATIONS,DUMMY,NDUMMY,CONVERGED,COORDS(:,MYNODE+1),DOSSTATS)
               NQ(MYNODE+1)=NQ(MYNODE+1)+1
               MINCOORDS(:,MYNODE+1)=COORDS(:,MYNODE+1) ! MINCOORDS contains quench coords
               VMINOLD=POTEL
               IF (CONVERGED.NE.1) WRITE(LFH, '(A)') 'bspt> WARNING - quench did not converge' 
               IF (DEBUG) WRITE(LFH, '(A,G20.10)') 'bspt> VMINOLD set to ',VMINOLD
            ENDIF
         ENDIF
!
! If either the PE (VNEW) or quench energy (VMINNEW) is out of range of the histogram
! then we reject the step and recount the previous one.
! Only recount if the quench energy is out of range for MINDENSITYT.
!
         IF (BSPT.AND.((I.GT.NEQUIL+PTSTEPS).OR.(MINDENSITYT)).AND.(MOD(I,1.0D0*QUENCHFRQ).EQ.0.0D0)) THEN
            BININDEX=INT((VMINNEW-HISTMIN)/HISTINT)+1
            IF ((BININDEX.GT.HBINMAX).OR.(BININDEX.LT.HBINMIN)) THEN
               RECOUNT=.TRUE.
               NOUTQBIN=NOUTQBIN+1
            ENDIF
         ENDIF
         IF (.NOT.MINDENSITYT) THEN
            IBININDEX=INT((VNEW-PTEMIN)/PEINT)+1
            IF ((IBININDEX.GT.NENRPER).OR.(IBININDEX.LT.1)) THEN
               RECOUNT=.TRUE.
               NOUTPEBIN=NOUTPEBIN+1
            ENDIF
         ELSE
            WCOMP=(VMINNEW-VMINOLD)*BETA(MYNODE) ! use difference in quench energies for MINDENSITYT
            W=MIN(1.0D0,EXP(-WCOMP))
            RANDOM=DPRAND()
            IF (RANDOM.GT.W) RECOUNT=.TRUE. 
         ENDIF

         VMINNEWSAVE=VMINNEW
         LBFGS_ITERATIONSSAVE=LBFGS_ITERATIONS
         VNEWSAVE=VNEW        ! this value is saved so it can be printed if the step is rejected
         IF ((CONVERGED.NE.1).OR.EVAPREJECT) RECOUNT=.TRUE. ! reject and recount
         IF (RECOUNT) THEN ! reject move
            DO K=1, NATOMS
               X(K)=XO(K)
               Y(K)=YO(K)
               Z(K)=ZO(K)
            ENDDO
            DO K=1,NATOMS
              COORDS(3*(K-1)+1,MYNODE+1)=X(K)
              COORDS(3*(K-1)+2,MYNODE+1)=Y(K)
              COORDS(3*(K-1)+3,MYNODE+1)=Z(K)
            ENDDO
!           WRITE(LFH,'(A,2G20.10)') 'bspt> recounting VOLD,VNEW=',VOLD,VNEW                
!           WRITE(LFH,'(A,2G20.10)') 'bspt> recounting VMINOLD,VMINNEW=',VMINOLD,VMINNEW
            VMINNEW=VMINOLD
            VNEW=VOLD
            DOSSTATS(1:MAXIT,1:2)=DOSSTATSO(1:MAXIT,1:2)
            LBFGS_ITERATIONS=LBFGS_ITERATIONSO
            IF (.NOT.FITS) NOUT(MYNODE)=NOUT(MYNODE)+1
         ELSE ! accept move
            NACCEPTPT(MYNODE)=NACCEPTPT(MYNODE)+1.0D0
            IACCEPT(MYNODE)=IACCEPT(MYNODE)+1
         ENDIF
         IF (.NOT.RECOUNT) THEN
            IF (MOD(I-1.0D0*1,1.0D0*PRTFRQ).EQ.1.0D0*0) WRITE(LFH, '(F15.1,A,G16.6,A,G16.6,A,G16.6,A,G16.6,A,I6,A)') 
     &        I,' En= ', VMINNEWSAVE, ' Eo= ',VMINOLD, ' Vn=', VNEWSAVE,' Vo=',VOLD,' iter ',LBFGS_ITERATIONSSAVE,' ACC'
         ELSE
            IF (MOD(I-1.0D0*1,1.0D0*PRTFRQ).EQ.1.0D0*0) WRITE(LFH, '(F15.1,A,G16.6,A,G16.6,A,G16.6,A,G16.6,A,I6,A)') 
     &        I,' En= ', VMINNEWSAVE, ' Eo= ',VMINOLD, ' Vn=', VNEWSAVE,' Vo=',VOLD,' iter ',LBFGS_ITERATIONSSAVE,' REJ'
         ENDIF
         IF (DEBUG) THEN
            CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
            WRITE(LFH,'(A,G20.10)') 'bspt> potential energy after acc/rej=',POTEL
         ENDIF
         CALL FLUSH(LFH)
!
!  At this point the quench and instantaneous energies for the current
!  configuration in the Markov chain are VMINNEW and VNEW.
!
!  QUENCHFRQ must be initialised to one for PTMC to avoid division by zero!
!
         IF (MOD(I,1.0D0*QUENCHFRQ).EQ.1.0D0*0) BININDEX=INT((VMINNEW-HISTMIN)/HISTINT)+1
         IBININDEX=INT((VNEW-PTEMIN)/PEINT)+1
!
!  Must not accumulate statistics until we have equilibrated for NEQUIL steps.
!
         IF (I.GT.NEQUIL) THEN
            IF (BSPT.AND.(MOD(I,1.0D0*QUENCHFRQ).EQ.0.0D0).AND.(I.GT.NEQUIL+PTSTEPS)) THEN 
               IF (BINSTRUCTURES.AND.(MOD(I, 1.0D0*SAVENTH).EQ.1.0D0*0).AND.(.NOT.RECOUNT)) THEN 
                  CALL SAVEBINSTRUCTURESMPI(POTEL,MINCOORDS(:,MYNODE+1),BININDEX,.TRUE.,MYNODE,NEWENERGY)
!                 IF (NEWENERGY) THEN 
!                     MINIMANUMBER(BININDEX, MYNODE)=MINIMANUMBER(BININDEX, MYNODE)+1
!                 ENDIF
               ENDIF
            ENDIF
            IF (BSPT.AND.(MOD(I,1.0D0*QUENCHFRQ).EQ.0.0D0).AND.(I.GT.NEQUIL+PTSTEPS)) THEN
               IF ((BININDEX.GT.HBINS).OR.(BININDEX.LT.1)) GOTO 888
            ENDIF
            IF ((IBININDEX.GT.NENRPER).OR.(IBININDEX.LT.1)) GOTO 888

            PEVISITS(IBININDEX,MYNODE)=PEVISITS(IBININDEX,MYNODE)+1.0D0
            IF (BSPT.AND.(MOD(I,1.0D0*QUENCHFRQ).EQ.0.0D0).AND.(I.GT.NEQUIL+PTSTEPS)) THEN ! otherwise BININDEX could be out of range
               QVISITS(BININDEX,MYNODE)=QVISITS(BININDEX,MYNODE)+1.0D0
               PEVISITS2(IBININDEX,BININDEX,MYNODE)=PEVISITS2(IBININDEX,BININDEX,MYNODE)+1
            ENDIF
            IF (BSPT.AND.(MOD(I,1.0D0*QUENCHFRQ).EQ.1.0D0*0).AND.(I.GT.NEQUIL+PTSTEPS)) THEN
               DISTANCEOLD=CALCULATEDDISTANCE(COORDS(:,MYNODE+1), COORDSO(:,MYNODE+1))
               DISTANCE(BININDEX,MYNODE)=DISTANCE(BININDEX,MYNODE)+DISTANCEOLD
            ENDIF
         ENDIF
!
! If we are not calculating distances then the following lines may be unnecessary.
!
888      DO K=1,NATOMS ! COORDS have to be resaved because we have to use noreset to be able to calc distance 
            COORDS(3*(K-1)+1,MYNODE+1)=X(K)
            COORDS(3*(K-1)+2,MYNODE+1)=Y(K)
            COORDS(3*(K-1)+3,MYNODE+1)=Z(K)
         ENDDO
!
!  VMINOLD is now set to VMINNEW and VOLD to VNEW.
!
         VOLD=VNEW
         VMINOLD=VMINNEW 

         DOSSTATSO(1:MAXIT,1:2)=DOSSTATS(1:MAXIT,1:2)
         LBFGS_ITERATIONSO=LBFGS_ITERATIONS

         IF ((I.LE.NEQUIL).AND.(MOD(I,1.0D0*NUPDATE).EQ.1.0D0*0)) THEN ! update MC step size if not fixed
            WAC=1.0*IACCEPT(MYNODE)/NUPDATE
            IF (.NOT.(FIXSTEP(MYNODE+1).OR.FIXBOTH(MYNODE+1))) THEN
!              IF (WAC.LT.0.4) THEN
               IF (WAC.LT.ACCRAT(MYNODE+1)-0.1D0) THEN
                  STEP(MYNODE+1)=STEP(MYNODE+1)*0.9D0
               ENDIF
!              IF (WAC.GT.0.6) THEN
               IF (WAC.GT.ACCRAT(MYNODE+1)+0.1D0) THEN
                  STEP(MYNODE+1)=STEP(MYNODE+1)*1.1D0
               ENDIF
            ENDIF
            IACCEPT(MYNODE)=0
            WRITE(LFH,'(A,G20.10)') 'bspt> maximum step size is now ',STEP(MYNODE+1)
         ENDIF
!
!  Replica exchange part follows. 

         E=VNEW
         IF (MYNODE.EQ.0) THEN
            RANDOM=DPRAND()
            J=(NPAR-1)*RANDOM
            RANDOM=DPRAND()
            IF (RANDOM.GT.EXCHPROB) THEN 
               J=-2
               EXCHANGE=.FALSE.
            ELSE
               EXCHANGE=.TRUE.
            ENDIF
         ENDIF
         CALL MPI_BCAST(J,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPIERR)
!
! J is set to -2 for no exchange. Should really just skip the next part?!
! We try to exchange with trajectory J+1. J starts from 0.
!
         IF (MYNODE.EQ.(J+1)) THEN
            CALL MPI_SEND(E,1,MPI_DOUBLE_PRECISION,J,0,MPI_COMM_WORLD,MPIERR)
         ENDIF
         IF (MYNODE.EQ.J) THEN
            CALL MPI_RECV(ER,1,MPI_DOUBLE_PRECISION,J+1,0,MPI_COMM_WORLD,IS,MPIERR)
            DBETA=BETA(J)-BETA(J+1)
            DELTA=E-ER
            W=MIN(1.0D0,DEXP(DELTA*DBETA))
            NTOT=NTOT+1 
            RANDOM=DPRAND()
            LESAVE=E
            IF (W.GT.RANDOM) THEN
               EXCHANGEACCEPT=.TRUE.
               IMESG=1
               CALL MPI_SEND(IMESG,1,MPI_INTEGER,J+1,0,MPI_COMM_WORLD,MPIERR)
               CALL MPI_SEND(ITRAJ,1,MPI_INTEGER,J+1,0,MPI_COMM_WORLD,MPIERR)
               CALL MPI_SEND(X,NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,MPIERR)
               CALL MPI_SEND(Y,NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,MPIERR)
               CALL MPI_SEND(Z,NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,MPIERR)
               CALL MPI_SEND(VMINNEW,1,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,MPIERR)
               CALL MPI_SEND(VNEW,1,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,MPIERR)
               CALL MPI_SEND(DOSSTATS,2*MAXIT,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,MPIERR)
               CALL MPI_SEND(LBFGS_ITERATIONS,1,MPI_INTEGER,J+1,1,MPI_COMM_WORLD,MPIERR)
               CALL MPI_RECV(ITRAJ,1,MPI_INTEGER,J+1,0,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(X,NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(Y,NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(Z,NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(VMINNEW,1,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(VNEW,1,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(DOSSTATS,2*MAXIT,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(LBFGS_ITERATIONS,1,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,IS,MPIERR)
               E=ER
               NEACCEPT=NEACCEPT+1
               DOSSTATSO(1:MAXIT,1:2)=DOSSTATS(1:MAXIT,1:2)
               LBFGS_ITERATIONSO=LBFGS_ITERATIONS
               VMINOLD=VMINNEW ! these variables are the same at this point in the program
               VOLD=VNEW
            ELSE
               EXCHANGEACCEPT=.FALSE.
               IMESG=0
               CALL MPI_SEND(IMESG,1,MPI_INTEGER,J+1,0,MPI_COMM_WORLD,MPIERR)
            ENDIF
            IF (EXCHANGEACCEPT) THEN
               IF (MOD(I-1.0D0,1.0D0*PRTFRQ).EQ.0) WRITE(LFH,'(A,I5,A,G16.6,A,I5,A,G16.6,A)') 
     &                        'bspt> Exchange from this replica ',J,' energy ',LESAVE,
     &                                                         ' to replica ',J+1,' energy ',ER,' ACC'
               IF (DEBUG) THEN
                  WRITE(LFH,'(A,G20.10)') 'bspt> Instantaneous PE is now ',VOLD

                  DO K=1,NATOMS
                     COORDS(3*(K-1)+1,MYNODE+1)=X(K)
                     COORDS(3*(K-1)+2,MYNODE+1)=Y(K)
                     COORDS(3*(K-1)+3,MYNODE+1)=Z(K)
                  ENDDO
                  CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
                  WRITE(LFH,'(A,G20.10)') 'bspt> Recalculated PE is ',POTEL,' VNEW=',VNEW
               ENDIF

            ELSE
               IF (MOD(I-1.0D0,1.0D0*PRTFRQ).EQ.0.0D0) WRITE(LFH,'(A,I5,A,G16.6,A,I5,A,G16.6,A)') 
     &                    'bspt> Exchange from this replica ',J,' energy ',E,' to replica ',
     &                           J+1,' energy ',ER,' REJ'
            ENDIF
         ENDIF
         CALL FLUSH(LFH)
         IF (MYNODE.EQ.(J+1)) THEN
            CALL MPI_RECV(IMESG,1,MPI_INTEGER,J,0,MPI_COMM_WORLD,IS,MPIERR)
            NTOT=NTOT+1
            ER=E
            IF (IMESG.EQ.1) THEN
!
!  Here we receive first, so we need to save variables that would otherwise be overwritten
!  before they are sent to the other replica. Or use XO, YO, ZO in the case of X, Y, Z.
!
               VOLDSAVE=VOLD
               DOSSTATSO(1:MAXIT,1:2)=DOSSTATS(1:MAXIT,1:2)
               LBFGS_ITERATIONSO=LBFGS_ITERATIONS
               CALL MPI_RECV(ITRAJO,1,MPI_INTEGER,J,0,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(XO,NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(YO,NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(ZO,NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(VMINOLD,1,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(VOLD,1,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(DOSSTATS,2*MAXIT,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_RECV(LBFGS_ITERATIONS,1,MPI_INTEGER,J,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_SEND(ITRAJ,1,MPI_INTEGER,J,0,MPI_COMM_WORLD,MPIERR)
               CALL MPI_SEND(X,NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,MPIERR)
               CALL MPI_SEND(Y,NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,MPIERR)
               CALL MPI_SEND(Z,NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,MPIERR)
               CALL MPI_SEND(VMINNEW,1,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_SEND(VNEW,1,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_SEND(DOSSTATSO,2*MAXIT,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
               CALL MPI_SEND(LBFGS_ITERATIONSO,1,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
               DO K=1, NATOMS
                  X(K)=XO(K)
                  Y(K)=YO(K)
                  Z(K)=ZO(K)
               ENDDO
               ITRAJ=ITRAJO
               NEACCEPT=NEACCEPT+1
               DOSSTATSO(1:MAXIT,1:2)=DOSSTATS(1:MAXIT,1:2)
               LBFGS_ITERATIONSO=LBFGS_ITERATIONS
               IF (MOD(I-1.0D0,1.0D0*PRTFRQ).EQ.0) WRITE(LFH,'(A,I5,A,G16.6,A,I5,A,G16.6,A)') 
     &                     'bspt> Exchange from replica ',MYNODE-1,' energy ',VOLD,
     &                                                         ' to this replica ',MYNODE,' energy ',VNEW,' ACC'
               VMINNEW=VMINOLD ! these variables are the same at this point in the program
               VNEW=VOLD
               IF (DEBUG) THEN
                  WRITE(LFH,'(2(A,G20.10))') 'bspt> Instantaneous PE is now ',VOLD,' and VNEW=',VNEW

                  DO K=1,NATOMS
                     COORDS(3*(K-1)+1,MYNODE+1)=X(K)
                     COORDS(3*(K-1)+2,MYNODE+1)=Y(K)
                     COORDS(3*(K-1)+3,MYNODE+1)=Z(K)
                  ENDDO
                  CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
                  WRITE(LFH,'(A,G20.10)') 'bspt> Recalculated PE is ',POTEL
               ENDIF

            ELSE
               IF (MOD(I-1.0D0,1.0D0*PRTFRQ).EQ.0.0D0) WRITE(LFH,'(A,I5,A,I5,A,G16.6,A)') 
     &                   'bspt> Exchange from      replica ',MYNODE-1,
     &                   ' energy not available to this replica ',MYNODE,' energy ',VOLD,' REJ'
            ENDIF
         ENDIF

         IF (I.GT.NEQUIL) THEN 
            EAV(MYNODE)=EAV(MYNODE)+E
            EAV2(MYNODE)=EAV2(MYNODE)+E**2
            IQE=INT((E-PTEMIN)/DHISTE+1)
            IF (IQE.GT.0.AND.IQE.LT.NHISTE) THEN
               NHISTQE(IQE,MYNODE)=NHISTQE(IQE,MYNODE)+1
            ENDIF
            DO K=1,NATOMS
               Q(1,K)=X(K)
               Q(2,K)=Y(K)
               Q(3,K)=Z(K)
            ENDDO
            IF (PERIODIC) THEN
               CALL QORDER_BLJ(Q,Q4,Q6)
            ELSE IF (CHRMMT) THEN
                 IF (ODIHET) CALL CHCALCDIHE(DIHEORDERPARAM,COORDS(1:3*NATOMS,MYNODE+1))
                 IF (OSASAT) CALL ORDER_SASA(SASAORDERPARAM,RPRO,COORDS(1:3*NATOMS:3,MYNODE+1),
     1                   COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1))
                 IF (ORGYT) CALL CHCALCRGYR(RGYR,COORDS(1:3*NATOMS:3,MYNODE+1),
     1                   COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1),.FALSE.) 
                 IF (OEINTT) CALL CHCALCEINT(EINT,COORDS(1:3*NATOMS:3,MYNODE+1),
     1                   COORDS(2:3*NATOMS:3,MYNODE+1),COORDS(3:3*NATOMS:3,MYNODE+1),POTEL)
            ELSE
               CALL QORDER_LJ(Q,Q4,Q6)
            ENDIF
            Q4AV(MYNODE)=Q4AV(MYNODE)+Q4
            Q4AV2(MYNODE)=Q4AV2(MYNODE)+Q4**2
            Q6AV(MYNODE)=Q6AV(MYNODE)+Q6
            Q6AV2(MYNODE)=Q6AV2(MYNODE)+Q6**2

            DIHEORDERPARAM_AV(MYNODE)=DIHEORDERPARAM_AV(MYNODE)+DIHEORDERPARAM
            DIHEORDERPARAM_AV2(MYNODE)=DIHEORDERPARAM_AV2(MYNODE)+DIHEORDERPARAM**2
            SASAORDERPARAM_AV(MYNODE)=SASAORDERPARAM_AV(MYNODE)+SASAORDERPARAM
            SASAORDERPARAM_AV2(MYNODE)=SASAORDERPARAM_AV2(MYNODE)+SASAORDERPARAM**2
            RGYR_AV(MYNODE)=RGYR_AV(MYNODE)+RGYR
            RGYR_AV2(MYNODE)=RGYR_AV2(MYNODE)+RGYR**2
            EINT_AV(MYNODE)=EINT_AV(MYNODE)+EINT
            EINT_AV2(MYNODE)=EINT_AV2(MYNODE)+EINT**2
            IF (CHRMMT) THEN
               IF (ODIHET) THEN
                  IQ4=INT((DIHEORDERPARAM-DIHEORDERPARAM_MIN)/DDIHE+1)
                  IF (IQ4.GT.0.AND.IQ4.LT.NHIST) THEN
                     NHISTQ4(IQ4,MYNODE)=NHISTQ4(IQ4,MYNODE)+1
                  ENDIF
                  IQ6=INT((SASAORDERPARAM-SASAORDERPARAM_MIN)/DSASA+1)             
                  IF (IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
                     NHISTQ6(IQ6,MYNODE)=NHISTQ6(IQ6,MYNODE)+1
                  ENDIF
               ELSE IF (ORGYT) THEN
                  IQ4=INT((RGYR-RGYR_MIN)/DRGYR+1)
                  IF (IQ4.GT.0.AND.IQ4.LT.NHIST) THEN
                     NHISTQ4(IQ4,MYNODE)=NHISTQ4(IQ4,MYNODE)+1
                  ENDIF
                  IQ6=INT((EINT-EINT_MIN)/DEINT+1)             
                  IF (IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
                     NHISTQ6(IQ6,MYNODE)=NHISTQ6(IQ6,MYNODE)+1
                  ENDIF
               ENDIF 
            ELSE
               IQ4=INT(Q4/DQ4+1)
               IF (IQ4.GT.0.AND.IQ4.LT.NHIST) THEN
                  NHISTQ4(IQ4,MYNODE)=NHISTQ4(IQ4,MYNODE)+1
               ENDIF
               IQ6=INT(Q6/DQ6+1)             
               IF (IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
                  NHISTQ6(IQ6,MYNODE)=NHISTQ6(IQ6,MYNODE)+1
               ENDIF
            ENDIF

! Free energy statistic

            IF (IQ4.GT.0.AND.IQ4.LT.NHIST.AND.IQ6.GT.0.AND.IQ6.LT.NHIST) THEN
               NHISTALLQ(IQ4,IQ6,MYNODE)=NHISTALLQ(IQ4,IQ6,MYNODE)+1
            ENDIF
            IF (IQ4.GT.0.AND.IQ4.LT.NHIST.AND.IQ6.GT.0.AND.IQ6.LT.NHIST.AND.
     1         IQE.GT.0.AND.IQE.LT.NHISTE) THEN
               NHISTGRAND(IQE,IQ4,IQ6,MYNODE)=NHISTGRAND(IQE,IQ4,IQ6,MYNODE)+1
            ENDIF
              
            IENR=IENR+1
!
! Dump visits histograms. 
! 
            IF (BSPTDUMPFRQ.GT.0) THEN
               IF (MOD(I,1.0D0*BSPTDUMPFRQ).EQ.0.0D0) THEN
                  WRITE (ISTR, '(I2,A1)') MYNODE, '.'
                  WRITE (SDUMMY, '(F15.1)') I
                  ISTR=TRIM(ADJUSTL(ISTR)) // TRIM(ADJUSTL(SDUMMY))
                  IF (I.GT.NEQUIL) THEN
                     FILENAME101="Visits.his."//TRIM(ADJUSTL(ISTR))
                     OPEN(UNIT=1986+MYNODE,FILE=FILENAME101, STATUS="unknown", form="formatted")
                     WRITE(1986+MYNODE, '(G20.10)') TEMPTRAJ(MYNODE)
                     WRITE(1986+MYNODE, '(A)') 'Visits to instantaneous PE bins without quench contributions'
                     DO K=1, NENRPER
                        WRITE(1986+MYNODE, '(G20.10,F20.1)') PTEMIN+(K-1)*PEINT,PEVISITS(K,MYNODE)
                     ENDDO
                     WRITE(1986+MYNODE, '(A)') 'Visits to quench bins'
                     DO K=1, HBINS
                          WRITE(1986+MYNODE, '(G20.10,F20.1)') BINLABEL(K), QVISITS(K,MYNODE)
                     ENDDO
                     CLOSE(1986+MYNODE)
                  ENDIF

                  IF (BSPT.AND.(I.GT.NEQUIL+PTSTEPS).AND.(.NOT.MINDENSITYT)) THEN
                     FILENAME101="Visits2.his."//TRIM(ADJUSTL(ISTR))
                     OPEN(UNIT=1986+MYNODE,FILE=FILENAME101,STATUS="unknown",FORM="UNFORMATTED")
                     WRITE(1986+MYNODE) PEVISITS2(1:NENRPER,1:HBINS,MYNODE)
                     CLOSE(1986+MYNODE)
                  ENDIF
!
!  Restart information
!
                  WRITE (ISTR,'(I2)') MYNODE
                  FILENAME9="bsptrestart."//TRIM(ADJUSTL(ISTR))
                  OPEN(UNIT=1985+MYNODE,FILE=FILENAME9, STATUS="unknown", form="formatted")
                  WRITE(1985+MYNODE,'(F20.1,3G20.10,F15.1,5I15)') I,VOLD,VMINOLD,STEP(MYNODE+1),
     &                                                NACCEPTPT(MYNODE),NEACCEPT,NTOT,NOUTQBIN,NOUTPEBIN,NQ(MYNODE+1)
                  DO J1=1,NATOMS
                     WRITE(1985+MYNODE,'(3G25.15)') X(J1),Y(J1),Z(J1)
                  ENDDO
                  CLOSE(1985+MYNODE)
               ENDIF
            ENDIF
         ENDIF
      ENDDO 
 
! End of main loop over BSPT or PT steps.
 
! compute the averages

      EAV(MYNODE)=EAV(MYNODE)/PTSTEPS
      EAV2(MYNODE)=EAV2(MYNODE)/PTSTEPS
      CV=(EAV2(MYNODE)-EAV(MYNODE)**2)*BETA(MYNODE)**2

      WRITE (ISTR, '(i10)') MYNODE+1
      FILENAME3="T.Ev.Cv.Ev2.Steps."//trim(adjustl(istr))
      OPEN(UNIT=41+MYNODE,FILE=FILENAME3, STATUS="unknown", form="formatted")
      IF (CHRMMT) THEN
         WRITE(41+MYNODE,'(6G20.10)') TEMPTRAJ(MYNODE),TEMPTRAJ(MYNODE)/0.001987, 
     &                                EAV(MYNODE), CV, EAV2(MYNODE),PTSTEPS
      ELSE
         WRITE(41+MYNODE,'(5G20.10)') TEMPTRAJ(MYNODE),EAV(MYNODE),CV, EAV2(MYNODE),PTSTEPS
      ENDIF
      CALL FLUSH(41+MYNODE)
      CLOSE(41+MYNODE)

      Q4AV(MYNODE)=Q4AV(MYNODE)/PTSTEPS
      Q6AV(MYNODE)=Q6AV(MYNODE)/PTSTEPS
      Q4AV2(MYNODE)=Q4AV2(MYNODE)/PTSTEPS
      Q6AV2(MYNODE)=Q6AV2(MYNODE)/PTSTEPS
      DIHEORDERPARAM_AV(MYNODE)=DIHEORDERPARAM_AV(MYNODE)/PTSTEPS
      SASAORDERPARAM_AV(MYNODE)=SASAORDERPARAM_AV(MYNODE)/PTSTEPS
      DIHEORDERPARAM_AV2(MYNODE)=DIHEORDERPARAM_AV2(MYNODE)/PTSTEPS
      SASAORDERPARAM_AV2(MYNODE)=SASAORDERPARAM_AV2(MYNODE)/PTSTEPS
      RGYR_AV(MYNODE)=RGYR_AV(MYNODE)/PTSTEPS
      EINT_AV(MYNODE)=EINT_AV(MYNODE)/PTSTEPS
      RGYR_AV2(MYNODE)=RGYR_AV2(MYNODE)/PTSTEPS
      EINT_AV2(MYNODE)=EINT_AV2(MYNODE)/PTSTEPS
      FQ4=SQRT(Q4AV2(MYNODE)-Q4AV(MYNODE)**2)
      FQ6=SQRT(Q6AV2(MYNODE)-Q6AV(MYNODE)**2)
      IF (BINARY) THEN  
         FILENAME4="T.Q4Av.Q6Av.Q4Av2.Q6Av2.Steps."//trim(adjustl(istr))
         OPEN(UNIT=1980,FILE=FILENAME4, STATUS="unknown", form="formatted")
         WRITE(1980,'(6G20.10)') TEMPTRAJ(MYNODE), Q4AV(MYNODE), Q6AV(MYNODE), Q4AV2(MYNODE),Q6AV2(MYNODE), PTSTEPS
         CLOSE(1980)
      ENDIF
      IF (CHRMMT) THEN 
         IF (ODIHET) THEN 
            FILENAME4="T.DiheAv.SasaAv.DiheAv2.SasaAv2.Steps."//trim(adjustl(istr))
            OPEN(UNIT=1980,FILE=FILENAME4, STATUS="unknown", form="formatted")
            WRITE(1980,'(6G20.10)') TEMPTRAJ(MYNODE), DIHEORDERPARAM_AV(MYNODE), 
     &              SASAORDERPARAM_AV(MYNODE), DIHEORDERPARAM_AV2(MYNODE),SASAORDERPARAM_AV2(MYNODE), PTSTEPS
            CLOSE(1980)
            ELSE IF (ORGYT) THEN
            FILENAME4="T.RgyrAv.EintAv.RgyrAv2.Eint2.Steps."//trim(adjustl(istr))
            OPEN(UNIT=1980,FILE=FILENAME4, STATUS="unknown", form="formatted")
            WRITE(1980,'(6G20.10)') TEMPTRAJ(MYNODE), RGYR_AV(MYNODE),EINT_AV(MYNODE), RGYR_AV2(MYNODE),  
     &              EINT_AV2(MYNODE), PTSTEPS
            CLOSE(1980)
         ENDIF
      ENDIF

      FILENAME5="profile_E."//trim(adjustl(istr)) 
      OPEN(UNIT=1981,FILE=FILENAME5, STATUS="unknown", form="formatted")
      DO K=1,NHISTE
         WRITE(1981,'(2G20.10)') PTEMIN+(K-1)*DHISTE,NHISTQE(K,MYNODE)
      ENDDO
      CLOSE(1981)

      IF (BINARY) THEN
         FILENAME6="profile_Q4."//trim(adjustl(istr))
         FILENAME7="profile_Q6."//trim(adjustl(istr))
         OPEN(UNIT=1982,FILE=FILENAME6, STATUS="unknown", form="formatted")
         OPEN(UNIT=1983,FILE=FILENAME7, STATUS="unknown", form="formatted")
         DO K=1,NHIST
            WRITE(1982,'(2G20.10)') (K-1)*DQ4,NHISTQ4(K,MYNODE)
            WRITE(1983,'(2G20.10)') (K-1)*DQ6,NHISTQ6(K,MYNODE)
         ENDDO
         CLOSE(1982)
         CLOSE(1983)
      ELSE IF (CHRMMT) THEN
         IF (ODIHET) THEN
            FILENAME6="profile_dihe."//trim(adjustl(istr))
            FILENAME7="profile_sasa."//trim(adjustl(istr))
            OPEN(UNIT=1982,FILE=FILENAME6, STATUS="unknown", form="formatted")
            OPEN(UNIT=1983,FILE=FILENAME7, STATUS="unknown", form="formatted")
            DO K=1,NHIST
               WRITE(1982,'(2G20.10)') diheorderparam_min+(K-1)*Ddihe,NHISTQ4(K,MYNODE)
               WRITE(1983,'(2G20.10)') SASAorderparam_min+(K-1)*Dsasa,NHISTQ6(K,MYNODE)
            ENDDO
            CLOSE(1982)
            CLOSE(1983)
         ELSE IF (ORGYT) THEN
            FILENAME6="profile_rgyr."//trim(adjustl(istr))
            FILENAME7="profile_eint."//trim(adjustl(istr))
            OPEN(UNIT=1982,FILE=FILENAME6, STATUS="unknown", form="formatted")
            OPEN(UNIT=1983,FILE=FILENAME7, STATUS="unknown", form="formatted")
            DO K=1,NHIST
               WRITE(1982,'(2G20.10)') rgyr_min+(K-1)*Drgyr,NHISTQ4(K,MYNODE)
               WRITE(1983,'(2G20.10)') eint_min+(K-1)*Deint,NHISTQ6(K,MYNODE)
            ENDDO
            CLOSE(1982)
            CLOSE(1983)
         ENDIF
      ENDIF

      WRITE(LFH, '(A,F15.1,A,G15.1,A,F15.5,A)') 'bspt> ',NACCEPTPT(MYNODE), ' steps accepted out of ', 
     &            PTSTEPS+NEQUIL+NQUENCH, ' i.e. ',NACCEPTPT(MYNODE)*100.0D0/(PTSTEPS+NEQUIL+NQUENCH),'%'
      WRITE(LFH, '(A,G20.10)') 'bspt> Final stepsize ', STEP(MYNODE+1)
      WRITE(LFH, '(A,G20.10,A,G20.10)') 'bspt> ',NEACCEPT, ' PT exchanges accepted out of ', NTOT
      IF (BSPT) WRITE(LFH, '(A,G20.10,A,G20.10,A)') 'bspt> ',NOUTQBIN, ' quenches outside quench bin range: ',
     &                                      NOUTQBIN*100.0D0/(NQUENCH),' %'
      WRITE(LFH, '(A,I8,A,G20.10,A)') 'bspt> ',NOUTPEBIN, ' potential energies outside bin range: ',
     &                                      NOUTPEBIN*100.0D0/(PTSTEPS+NEQUIL+NQUENCH),' %'
      WRITE(LFH, '(A,I8)') 'bspt> total number of quenches actually required=',NQ(MYNODE+1)

      DO K=1,NATOMS
         Q(1,K)=X(K)
         Q(2,K)=Y(K)
         Q(3,K)=Z(K)
      ENDDO
      IF (PERIODIC) THEN
         CALL QORDER_BLJ(Q,Q4,Q6)
      ELSE
         CALL QORDER_LJ(Q,Q4,Q6)
      ENDIF

!
! Dump visits histograms. Energies are written for the middle of the bins.
! This is the really important part!
!
      WRITE (ISTR, '(i10)') MYNODE 
      FILENAME101="Visits.his."//trim(adjustl(istr))
      OPEN(UNIT=1986+MYNODE,FILE=FILENAME101, STATUS="unknown", form="formatted")
      WRITE(1986+MYNODE, '(G20.10)') TEMPTRAJ(MYNODE)
      WRITE(1986+MYNODE, '(A)') 'Visits to instantaneous PE bins without quench contributions'
      DO K=1, NENRPER
         WRITE(1986+MYNODE, '(G20.10,F20.1)') PTEMIN+(K-1)*PEINT,PEVISITS(K,MYNODE)
      ENDDO
      WRITE(1986+MYNODE, '(A)') 'Visits to quench bins'
      DO K=1, HBINS
           WRITE(1986+MYNODE, '(G20.10,F20.1)') BINLABEL(K), QVISITS(K,MYNODE)
      ENDDO
      CLOSE(1986+MYNODE)

      IF (BSPT.AND.(.NOT.MINDENSITYT)) THEN
         FILENAME101="Visits2.his."//trim(adjustl(istr))
         OPEN(UNIT=1986+MYNODE,FILE=FILENAME101,STATUS="unknown",FORM="UNFORMATTED")
         WRITE(1986+MYNODE) PEVISITS2(1:NENRPER,1:HBINS,MYNODE)
         CLOSE(1986+MYNODE)
      ENDIF
!
! Dump final restart information.
!
      WRITE (ISTR,'(I2)') MYNODE
      FILENAME9="bsptrestart."//TRIM(ADJUSTL(ISTR))
      OPEN(UNIT=1985+MYNODE,FILE=FILENAME9, STATUS="unknown", form="formatted")
      WRITE(1985+MYNODE,'(F20.1,3G20.10,F15.1,5I10)') I,VOLD,VMINOLD,STEP(MYNODE+1),
     &                                          NACCEPTPT(MYNODE),NEACCEPT,NTOT,NOUTQBIN,NOUTPEBIN,NQ(MYNODE+1)
      DO J1=1,NATOMS
         WRITE(1985+MYNODE,'(3G25.15)') X(J1),Y(J1),Z(J1)
      ENDDO
      CLOSE(1985+MYNODE)

      RETURN
!op226 End MPI section }}}
#else
      RETURN
#endif
      END SUBROUTINE PTBASINSAMPLING
