!  NEWRES {{{ 
!  Reseed if the energy has not improved by more than ECONV over the
!  last NRELAX mc steps.
!  If AVOID is true then save the energy and coordinates of the lowest
!  minimum achieved before each restart and restart if we get too close
!  to any one of them using bipartite matching and mind. Note that bipartite
!  matching and mind can give a local minimum of distance if the optimal 
!  permutation-inversion isn;t found. Using ORIENT to effect a standard
!  orientation first seems to help. It should at least ensure that permutation-inversion
!  isomers are always found. 
!  Would it be possible to just use the energy and inertia components to identify
!  permutation-inversion isomers instead of MINPERM and MIND?
!  This would be like using a small value of AVOIDDIST, which doesn;t seem to be
!  as good.
!
      SUBROUTINE NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,&
     &                  RMIN,RVAT,BRUN,SCREENC,QDONE,JACCPREV,NSUCCESS,NFAIL,NFAILT,NSUCCESST)
      USE commons
      !USE modamber9, only : mdstept 
      
      IMPLICIT NONE
      INTEGER J1, JP, JBEST(NPAR), ITERATIONS, J2, JACCPREV, BRUN, QDONE, J3, PERM(NATOMS), NPERM, NTRIES
      INTEGER, PARAMETER :: MAXIMUMTRIES=20
      DOUBLE PRECISION EBEST(NPAR), BESTCOORDS(3*NATOMS,NPAR)
      DOUBLE PRECISION :: EPPREV(NPAR), POTEL, TIME, RCOORDS(3*NATOMS), DIST2, DUMMY2
      DOUBLE PRECISION :: RVAT(NATOMS), RMIN, RANDOM, SR3, SCREENC(3*NATOMS)
      DOUBLE PRECISION :: DPRAND, FCOORDS(3*NATOMS), XMSBSAVE(3*NATOMS)
      DOUBLE PRECISION :: DUMMY(3*NATOMS), DISTANCE, XMSB(3*NATOMS)
      DOUBLE PRECISION :: EBESTP, BESTCOORDSP(3*NATOMS), WORSTRAD, RMAT(3,3), QENERGY
      DOUBLE PRECISION :: ROTA(3,3), ROTINVA(3,3)
      INTEGER NSUCCESS(NPAR), NFAIL(NPAR), NFAILT(NPAR), NSUCCESST(NPAR), NORBIT1, NORBIT2, INVERT, NPOSITION
      LOGICAL RES1, RES2, HIGHEST
      COMMON /MYPOT/ QENERGY

      SR3=DSQRT(3.0D0)
      IF (POTEL.LT.EBEST(JP)) THEN
         IF (EBEST(JP)-POTEL.GT.ECONV) JBEST(JP)=J1
         EBESTP=EBEST(JP)
         BESTCOORDSP(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP) ! save previous BESTCOORDS for possible use
         EBEST(JP)=POTEL ! reset ebest, but not necessarily jbest
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
      ENDIF
!
!  Reseed if the energy has not improved in the last NRELAX mc cycles,
!  or if the current minimum is too close to one of the NMSBSAVE structures
!  saved in MSBCOORDS.
!
!  Instead of using the current minimum, employ the current minimum in 
!  BESTCOORDS for the AVOID check. Then we need only do the AVOID check when we have
!  a new minimum in BESTCOORDS, i.e. J1.EQ.JBEST(JP).
!
      RES1=.FALSE.
      IF (J1-JBEST(JP).GT.NRELAX) RES1=.TRUE.
!     WRITE(LFH,'(A,I5,2G17.7,3I5,L8)') 'J1,POTEL,EBEST,JBEST,J1-JBEST,NRELAX,RES1=',
!    1                                J1,POTEL,EBEST(JP),JBEST(JP),J1-JBEST(JP),NRELAX,RES1
      RES2=.FALSE.
!     IF ((.NOT.RES1).AND.AVOID) THEN
!     PRINT*,'RES1,AVOID,J1,JBEST(JP)=',RES1,AVOID,J1,JBEST(JP)
      IF ((.NOT.RES1).AND.AVOID.AND.(J1.EQ.JBEST(JP)).AND.(NMSBSAVE.GT.0)) THEN ! best minimum has just changed.
         FCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
         savedloop: DO J2=1,MIN(NMSBSAVE,MAXSAVE)
!
!  Bipartite matching routine for permutations. Coordinates in FCOORDS do not change
!  but the coordinates in XMSB do. DISTANCE is the distance in this case.
!
! If the energy is lower no reseeding regardless of separation ? DJW
!           PRINT '(A,2G20.10,L8)','POTEL,MSBE(J2)-ECONV,POTEL.LT.MSBE(J2)-ECONV=',
!    1                              POTEL,MSBE(J2)-ECONV,POTEL.LT.MSBE(J2)-ECONV
            IF (POTEL.LT.MSBE(J2)-ECONV) CYCLE savedloop
            XMSB(1:3*NATOMS)=MSBCOORDS(1:3*NATOMS,J2)
            CALL MINPERMDIST(FCOORDS,XMSB,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DISTANCE,DIST2,RIGID,RMAT)
            IF (DISTANCE.LT.AVOIDDIST) THEN
               RES2=.TRUE.
               WRITE(LFH,'(A,G20.10,A,I6,A,G20.10,A,F10.3)') 'newres> Minimum energy ',POTEL,&
     &                                    ' is too close to saved structure ',&
     &                                 J2,' with energy ',MSBE(J2),' dist=',DISTANCE
               GOTO 20
            ENDIF
         ENDDO savedloop
      ENDIF

20    CONTINUE 


      IF (RES1.OR.RES2) THEN
          !  Does not seem necessary to save MSB data for restarts.
          !  If we are reseeding because RES2 is true then:
          !  (1) If the AVOID test is based on the distance for BESTCOORDS, add
          !      these to the  MSB list if the distance is > 0.01D0, i.e. a different structure.
          !      Otherwise, add the coordinates of the previous lowest energy minimum to
          !      the MSB list; these have been saved in BESTCOORDSP and the energy in EBESTP.
          !  (2) Alternatively, if the AVOID test is based on the coordinates of the current
          !      structure, we should add these to the MSB list instead. Testing
          !      every step involves an overhead, which grows with the number of saved taboo
          !      structures. We could make the taboo list cyclic, so that new structures to be
          !      avoided replace old ones at the beginning of the list when the length of the
          !      list is exceeded. We are currently only doing the bipartite matching test if
          !      the energy of the lowest minimum since the last reseeding has changed, which
          !      saves time.
          !
          !  Change to a cyclic AVOID list 30/12/04.

         IF (RES2.AND.(.NOT.AVOIDRESEEDT)) THEN 
            POTEL=MAX(1.0D10,10.0D0*POTEL) ! This should be enough to reject the step until we do a massive Thomson problem!
            WRITE(LFH,'(A,I8,A)') 'newres> Resetting energy to a large value due to taboo condition'
            RETURN
         ELSEIF (RES1.OR.(RES2.AND.(DISTANCE.GT.0.1D0))) THEN ! new condition
            NMSBSAVE=NMSBSAVE+1
            NPOSITION=MOD(NMSBSAVE,MAXSAVE)
            IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
            MSBE(NPOSITION)=EBEST(JP)
            FCOORDS(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP)

            CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG,ROTA,ROTINVA,STOCKT)

            MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)

            WRITE(LFH,'(A,I6,A,G20.10)') 'newres> Moving current best minimum to position ',NPOSITION,&
     &                         ' in the AVOID list E=',EBEST(JP)
!           OPEN(UNIT=34,FILE='MSBdata',POSITION='APPEND')
!           WRITE(34,'(G20.10)') MSBE(NPOSITION) 
!           WRITE(34,'(3G20.10)') BESTCOORDS(1:3*NATOMS,JP)
!           CLOSE(34)
         ELSEIF (RES2.AND.(DISTANCE.LE.0.01D0)) THEN
!           IF (NMSBSAVE.LT.MAXSAVE) THEN
               NMSBSAVE=NMSBSAVE+1
               NPOSITION=MOD(NMSBSAVE,MAXSAVE)
               IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
               MSBE(NPOSITION)=EBESTP
               FCOORDS(1:3*NATOMS)=BESTCOORDSP(1:3*NATOMS)
               CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG,ROTA,ROTINVA,STOCKT)
               MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)
               WRITE(LFH,'(A,I6,A,G20.10)') 'newres> Moving previous best minimum to position ',NPOSITION,&
     &                             ' in the AVOID list E=',EBESTP
!           ENDIF
         ENDIF ! end new condition

         IF (NEWRESTART_MD) THEN ! lb415
            IF (RES1) WRITE(LFH,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - perturbing'
            WRITE(LFH,'(A,I8,A)') 'newres> Reseeding via a short high temperature MD run'
            CHANGE_TEMP = .true.

! th368: 20-10-2009 Extending MD-Reseeding to the CHARMM interface
! terminate if neither AMBER nor CHARMM was requested in the data file

             IF (AMBERT) THEN
               CALL TAKESTEPAMBER(JP,COORDS(:,JP),movableatomlist,nmovableatoms,ligmovet,mdstept,randomseedt)
             ELSEIF (CHRMMT) THEN
               CALL CHMD(JP)
             ELSE
               WRITE(LFH,'(A,I8,A)') 'newres> Molecular Dynamics reseeding is available for AMBER or CHARMM runs only.'
               STOP
             ENDIF
! end th368: 20-10-2009

            CHANGE_TEMP = .false.
         ELSE
            IF (NHSRESTART.GT.0) THEN
               IF (RES1) WRITE(LFH,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - perturbing'
               IF (RES2) WRITE(LFH,'(A,I8,A)') 'newres> Reseeding due to taboo condition'
               IF (SHELLMOVES(JP)) WRITE(LFH,'(A)') 'newres> Turning off shell moves'
               SHELLMOVES(JP)=.FALSE.
               CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
            ELSE
               IF (RES1) THEN
                  IF (AVOIDRESEEDT) THEN
                     WRITE(LFH,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - reseeding'
                  ELSE
                     WRITE(LFH,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP)
                  ENDIF
               ENDIF
               IF (RES2) WRITE(LFH,'(A,I8,A)') 'newres> Reseeding due to taboo condition'
               HIGHEST=.TRUE.
               DO J2=1,NPAR
                  IF (J2.EQ.JP) CYCLE
                  IF (EBEST(J2).GT.EBEST(JP)) THEN
                     HIGHEST=.FALSE.
                     EXIT
                  ENDIF
               ENDDO
               IF (HIGHEST) THEN
                  WRITE(LFH,'(A,I6,A,F20.10)') 'newres> Parallel run ',JP,' has the highest energy ',EBEST(JP)
                  IF (NPAR.GT.1) WRITE(LFH,'(6F20.10)') EBEST(1:NPAR)
                  IF (.NOT.AVOIDRESEEDT) THEN
                     WRITE(LFH,'(A,I8,A)') 'newres> Resetting energy to a large value to reject this step (not reseeding)'
                     JBEST(JP)=J1
                     EBEST(JP)=POTEL ! this is communicated via common block MYPOT, which is in quench.f
                     POTEL=MAX(1.0D10,10.0D0*POTEL) ! This should be enough to reject the step until we do a massive Thomson problem!
                     RETURN
                  ELSE
                     WRITE(LFH,'(A)') 'newres> Full reseeding'
                     DO J2=1,3*NATOMS
                        RANDOM=(DPRAND()-0.5D0)*2.0D0
                        COORDS(J2,JP)=RANDOM*DSQRT(RADIUS)/SR3
                     ENDDO
                     NCORE(JP)=0
                     PTGROUP(JP)='   '
                  ENDIF
               ELSEIF (NCORE(JP).GT.0) THEN
!         
!  Reseed everything if this is the lowest of a set of parallel runs, otherwise
!  just reseed the surface.
!         
                  WRITE(LFH,'(A)') 'newres> Accepting an enforced surface reseeding'
                  DUMMY2=-1.0D0
                  DO J2=NATOMS-NCORE(JP)+1, NATOMS
                     DISTANCE=COORDS(3*(J2-1)+1,JP)**2+COORDS(3*(J2-1)+2,JP)**2+COORDS(3*(J2-1)+3,JP)**2
                     IF (DISTANCE.GT.DUMMY2) DUMMY2=DISTANCE
                  ENDDO
                  DISTANCE=SQRT(DISTANCE)
                  WRITE(LFH,'(A,F15.5)') 'newres> largest core radius=',DISTANCE
                  DO J2=1,NATOMS-NCORE(JP)
                     COORDS(3*(J2-1)+1,JP)=(DPRAND()-0.5D0)*2.0D0
                     COORDS(3*(J2-1)+2,JP)=(DPRAND()-0.5D0)*2.0D0
                     COORDS(3*(J2-1)+3,JP)=(DPRAND()-0.5D0)*2.0D0
                     RANDOM=DISTANCE+1.0D0+DPRAND()*(SQRT(RADIUS)-DISTANCE-1.0D0)
                     DUMMY2=SQRT(COORDS(3*(J2-1)+1,JP)**2+COORDS(3*(J2-1)+2,JP)**2+COORDS(3*(J2-1)+3,JP)**2)
                     COORDS(3*(J2-1)+1,JP)=COORDS(3*(J2-1)+1,JP)*RANDOM/DUMMY2
                     COORDS(3*(J2-1)+2,JP)=COORDS(3*(J2-1)+2,JP)*RANDOM/DUMMY2
                     COORDS(3*(J2-1)+3,JP)=COORDS(3*(J2-1)+3,JP)*RANDOM/DUMMY2
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
!
!  This step will be accepted because EPREV(JP)=POTEL. However, there
!  could be additional steps, so we need
!  to change COORDSO and VATO. Should reset EBEST and BESTCOORDS as well.
!
!  next line should be uncommented if routine is made available to use with CHARMM
!         IF (CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
654      CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         IF (RMS.GT.BQMAX) THEN
            WRITE(LFH,'(A)') 'newres> Quench from reseeded geometry failed - try again'
            DO J2=1,3*NATOMS
               RANDOM=(DPRAND()-0.5D0)*2.0D0
               COORDS(J2,JP)=RANDOM*DSQRT(RADIUS)/SR3
            ENDDO
            GOTO 654
         ENDIF
         NSUCCESS(JP)=0
         NFAIL(JP)=0
         POTEL=QENERGY
         EBEST(JP)=POTEL ! this is communicated via common block MYPOT, which is in quench.f
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         JBEST(JP)=J1
         EPREV(JP)=POTEL !
         EPPREV(JP)=0.0D0
         NSYMREM=0
! th368: 20.10.2009 updates should only take place in case
! reseeding took place - transfered from subroutine MC lines 363 to 367

         IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,JP))
         COORDSO(1:3*(NATOMS-NSEED),JP)=COORDS(1:3*(NATOMS-NSEED),JP)
!        WRITE(LFH,'(A,2G20.10)') 'newres> coordso changed: ',COORDSO(1,JP),COORDS(1,JP)
         VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)

! end th368 20.10.2009
      ENDIF

      RETURN
      END
      ! }}}
