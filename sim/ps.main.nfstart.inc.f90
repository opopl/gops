!
!  NFSTART is the first true frequency to use from a path.info file, and
!  NFFINISH is the last non-zero frequency to read for a minimum.
!  Used in getnewpath and tssearch. KAPPA is the number of non-zero 
!  vibrational frequencies.
!
NFSTART=1
      NFFINISH=3*NATOMS-6
      KAPPA=3*NATOMS-6
      NGLY=0
      IF (PULLT) THEN
         NFFINISH=3*NATOMS-4
         KAPPA=3*NATOMS-4
      ENDIF


