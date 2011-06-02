
WRITE(LFH,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
WRITE(LFH,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'

DO ISTEP=1,NSTEPS 
  MCTEMP = TEMP
  CALL TAKESTEP

  NQ=NQ+1
  CALL QUENCH(SCREENC,ITERATIONS,TIME,QDONE)  

  include mc.w.111.i.f90
  include mc.bh.pairdist.inc.f90
  include mc.bh.trackdata.inc.f90

  CALL TRAN(E,EPREV,ATEST,MCTEMP)

  include mc.bh.checkmarkov.inc.f90

  ! Accept or reject step
  include mc.bh.accrej.inc.f90
  ! Check the acceptance ratio every NACCEPT steps and update:
  !     NSUCCESST (the total number of accepted moves)
  !     NFAILT (the total number of rejected moves)

  IF ( (MOD(ISTEP,NACCEPT).EQ.0) CALL ACCREJ(NSUCCESS,NFAIL,NSUCCESST,NFAILT)
    TEMP=TEMP*SCALEFAC
  ENDIF
  IF (NQ.GT.NSTEPS) GOTO 37

ENDDO

