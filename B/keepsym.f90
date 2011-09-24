      SUBROUTINE KEEPSYM(NP)

      USE COMMONS

      IMPLICIT NONE
      INTEGER NP, J2, SYMINDEX(NATOMS), J3, J4
      DOUBLE PRECISION LCOORDS(3*NATOMS), NEWQ(3*NATOMS), SYMDELTA(3*NATOMS), DELTA(3*NATOMS),  SYMOP1(3,3)
      DOUBLE PRECISION STEPLENGTH, NEWSTEPLENGTH, ODIST1, ODIST2, DMIN, DUMMY
      LOGICAL ASSIGNED(NATOMS), BAD

      LCOORDS(1:3*NATOMS)=COORDSO(1:3*NATOMS,NP)
      DELTA(1:3*NATOMS)=COORDS(1:3*NATOMS,NP)-COORDSO(1:3*NATOMS,NP)
!     STEPLENGTH=SUM(DELTA(1:3*NATOMS)**2)
      SYMDELTA(1:3*NATOMS)=DELTA(1:3*NATOMS)

!     OPEN(UNIT=789,FILE='symstep.xyz',STATUS='UNKNOWN')
!     WRITE(789,'(I6)') NATOMS
!     WRITE(789,'(A)') 'coordso before symmetrising step'
!     WRITE(789,'(A3,3G20.10)') ('LA ',COORDSO(3*(J2-1)+1,NP),COORDSO(3*(J2-1)+2,NP),
!    &                                 COORDSO(3*(J2-1)+3,NP),J2=1,NATOMS-NCORE(NP))
!     WRITE(789,'(A3,3G20.10)') ('LB ',COORDSO(3*(J2-1)+1,NP),COORDSO(3*(J2-1)+2,NP),
!    &                                 COORDSO(3*(J2-1)+3,NP),J2=NATOMS-NCORE(NP)+1,NATOMS)
!     WRITE(789,'(I6)') NATOMS
!     WRITE(789,'(A)') 'coords before symmetrising step'
!     WRITE(789,'(A3,3G20.10)') ('LA ',COORDS(3*(J2-1)+1,NP),COORDS(3*(J2-1)+2,NP),
!    &                                 COORDS(3*(J2-1)+3,NP),J2=1,NATOMS-NCORE(NP))
!     WRITE(789,'(A3,3G20.10)') ('LA ',COORDS(3*(J2-1)+1,NP),COORDS(3*(J2-1)+2,NP),
!    &                                 COORDS(3*(J2-1)+3,NP),J2=NATOMS-NCORE(NP)+1,NATOMS)
!
!  New algorithm - choose the closest unclaimed atom in each case, so that
!  no tolerances are involved. 
!
      DO J2=1,NSYMREM
         BAD=.FALSE.
         SYMOP1(1:3,1:3)=SYMREM(J2,1:3,1:3)
         CALL MATMULV(NEWQ,LCOORDS,SYMOP1,NATOMS,3,3)
         ASSIGNED(1:NATOMS)=.FALSE.
         SYMINDEX(1:NATOMS)=0
         DO J3=1,NATOMS
            DMIN=1.0D100
            DO J4=1,NATOMS
               DUMMY=(LCOORDS(3*(J4-1)+1)-NEWQ(3*(J3-1)+1))**2+ &
     &               (LCOORDS(3*(J4-1)+2)-NEWQ(3*(J3-1)+2))**2+ &
     &               (LCOORDS(3*(J4-1)+3)-NEWQ(3*(J3-1)+3))**2
!              PRINT '(A,2I5,2G15.5)','J3,J4,DUMMY,DMIN=',J3,J4,DUMMY,DMIN
               IF (DUMMY.LT.DMIN) THEN
                  IF (ASSIGNED(J4)) THEN
!                    WRITE(LFH,'(2(A,I5),A,F12.3)') 'WARNING closest atom ',J4,' to image atom ',J3,
!    &                                       ' already assigned, dist=',SQRT(DUMMY)
                  ELSE
                     IF (SYMINDEX(J3).GT.0) ASSIGNED(SYMINDEX(J3))=.FALSE.
                     SYMINDEX(J3)=J4
                     ASSIGNED(J4)=.TRUE.
                     DMIN=DUMMY
                  ENDIF
               ENDIF
            ENDDO
            IF (DEBUG.AND.(SQRT(DMIN).GT.SYMTOL5)) THEN
               WRITE(LFH, '(2(A,I5),A,F12.3)') 'WARNING closest image to atom ',J3,' is ',SYMINDEX(J3), &
     &                                         ' distance=',SQRT(DMIN)
               BAD=.TRUE.
            ENDIF
         ENDDO
!        IF (BAD) THEN
!           OPEN(UNIT=1,FILE='keepsym.xyz',STATUS='UNKNOWN')
!           WRITE(1,*) NATOMS
!           WRITE(1,'(A)') 'LCOORDS'
!           DO J4=1,NATOMS
!              WRITE(1,'(A2,3X,3F20.10)') 'LA',LCOORDS(3*(J4-1)+1),LCOORDS(3*(J4-1)+2),LCOORDS(3*(J4-1)+3)
!           ENDDO
!           WRITE(1,*) NATOMS
!           WRITE(1,'(A)') 'NEWQ'
!           DO J4=1,NATOMS
!              WRITE(1,'(A2,3X,3F20.10)') 'LA',NEWQ(3*(J4-1)+1),NEWQ(3*(J4-1)+2),NEWQ(3*(J4-1)+3)
!           ENDDO
!           CLOSE(1)
!        ENDIF

         CALL MATMULV(NEWQ,DELTA,SYMOP1,NATOMS,3,3)
!
!  For each symmetry operation we average over the steps for the atoms
!  that are mapped on to after acting with the symmetry operation on the
!  original step vector.
!
         DO J3=1,NATOMS
            SYMDELTA(3*(J3-1)+1)=SYMDELTA(3*(J3-1)+1)+NEWQ(3*(SYMINDEX(J3)-1)+1)
            SYMDELTA(3*(J3-1)+2)=SYMDELTA(3*(J3-1)+2)+NEWQ(3*(SYMINDEX(J3)-1)+2)
            SYMDELTA(3*(J3-1)+3)=SYMDELTA(3*(J3-1)+3)+NEWQ(3*(SYMINDEX(J3)-1)+3)
         ENDDO
      ENDDO

!     NEWSTEPLENGTH=SUM(SYMDELTA(1:3*NATOMS)**2)
!     SYMDELTA(1:3*NATOMS)=SYMDELTA(1:3*NATOMS)*SQRT(STEPLENGTH/NEWSTEPLENGTH)
!
!  Maintain CofM distance in symmetry adapted step.
!
      SYMDELTA(1:3*NATOMS)=SYMDELTA(1:3*NATOMS)/(1+NSYMREM)

      DO J2=1,NATOMS
         ODIST1=SUM(COORDSO(3*(J2-1)+1:3*(J2-1)+3,NP)**2)
         ODIST2=SUM((COORDSO(3*(J2-1)+1:3*(J2-1)+3,NP)+SYMDELTA(3*(J2-1)+1:3*(J2-1)+3))**2)
         IF (ODIST2.NE.0.0D0) COORDS(3*(J2-1)+1:3*(J2-1)+3,NP)=(COORDSO(3*(J2-1)+1:3*(J2-1)+3,NP) &
     &        +SYMDELTA(3*(J2-1)+1:3*(J2-1)+3)) * SQRT(ODIST1/ODIST2)
      ENDDO

!     DO J2=1,NATOMS
!        DO J3=J2+1,NATOMS
!           DUMMY=SQRT((COORDS(3*(J3-1)+1,NP)-COORDS(3*(J2-1)+1,NP))**2+
!    &                 (COORDS(3*(J3-1)+2,NP)-COORDS(3*(J2-1)+2,NP))**2+
!    &                 (COORDS(3*(J3-1)+3,NP)-COORDS(3*(J2-1)+3,NP))**2)
!           IF (DUMMY.LT.0.1D0) WRITE(LFH,'(A,I6,A,I6,A,G20.10)') 'WARNING *** atoms ',J3,' and ',J2,
!    &             ' distance=',DUMMY
!        ENDDO
!     ENDDO

!     WRITE(789,'(I6)') NATOMS
!     WRITE(789,'(A)') 'after symmetrising step'
!     WRITE(789,'(A3,3G20.10)') ('LA ',COORDS(3*(J2-1)+1,NP),COORDS(3*(J2-1)+2,NP),
!    &                                 COORDS(3*(J2-1)+3,NP),J2=1,NATOMS-NCORE(NP))
!     WRITE(789,'(A3,3G20.10)') ('LA ',COORDS(3*(J2-1)+1,NP),COORDS(3*(J2-1)+2,NP),
!    &                                 COORDS(3*(J2-1)+3,NP),J2=NATOMS-NCORE(NP)+1,NATOMS)
!     CLOSE(789)
    
      RETURN

      END SUBROUTINE
