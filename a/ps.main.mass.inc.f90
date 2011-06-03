      MASS=1.0D0 ; ZSYMBOL=ZSYM
      IF (NTAG.GT.0) THEN
            MASS=MASS*TAGFAC
      ENDIF

!  If there is a file called mass in the current directory then use these masses instead
!  of the defaults.
      CALL INQF('mass',MASSFILE)
      IF (MASSFILE) THEN
         WRITE(*,'(A)') 'main> Reading individual atomic symbols and masses from file masses'
         CALL OPENF(1,'O','mass')
         READ (1,*) (ZSYMBOL(J1),MASS(J1),J1=1,NATOMS)
         CLOSE(1)
      ENDIF

