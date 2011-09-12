
WRITE(LFH,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
WRITE(LFH,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'

WRITE(*,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
WRITE(*,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'

! main basin hopping loop 

bhloop: DO ISTEP=NDONE+1,NSTEPS 

  MCTEMP = TEMP
  CALL TAKESTEP
  write(*,*) ISTEP,'  TAKESTEP done'

  NQ=NQ+1
  CALL QUENCH(SCREENC,ITERATIONS,TIME,CONVG)  
  write(*,*) ISTEP,'  QUENCH done'

  include "mc.w.111.i.f90"
  !include "mc.bh.pairdist.inc.f90"
  !include "mc.bh.trackdata.inc.f90"

  CALL TRAN(QE,QEPREV,ATEST,MCTEMP)

  !include "mc.bh.checkmarkov.inc.f90"

  ! Accept or reject step
  include "mc.bh.accrej.i.f90"
  ! Check the acceptance ratio every NACCEPT steps and update:
  !     NSUCCESST (the total number of accepted moves)
  !     NFAILT (the total number of rejected moves)

  IF ( (MOD(ISTEP,NACCEPT).EQ.0)) CALL ACCREJ(NSUCCESS,NFAIL,NSUCCESST,NFAILT)
    TEMP=TEMP*SCALEFAC
  IF (NQ.GT.NSTEPS) GOTO 37

ENDDO bhloop


