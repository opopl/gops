
      SUBROUTINE TAKESTEP(NP)

      USE COMMONS
      USE V
      USE F

      IMPLICIT NONE
      
      ! sub
      INTEGER,INTENT(IN) ::    NP 
      ! loc
      DOUBLE PRECISION,DIMENSION(NR) ::   RND
      DOUBLE PRECISION ::   DMAX,VMAX,VMIN,localstep
      INTEGER J1,J2

       IF (ABS(COORDSO(1,NP)-COORDS(1,NP)).GT.1.0D-3) THEN
          WRITE(LFH,'(A,2G20.10)'),'takestep> WARNING - coordso will be changed: ',COORDSO(1,NP),COORDS(1,NP)
       ENDIF

      COORDSO=COORDS; VATO=VAT

        LOCALSTEP=STEP(NP)
        !CALL GETSEED(THE_SEED)
        !CALL GETRND(RND,NR,-1.0D0,1.0D0,THE_SEED)
        CALL GETRND(RND,NR,-1.0D0,1.0D0)
        COORDS(1:NR,NP)=COORDS(1:NR,NP)+LOCALSTEP*RND(1:NR)

      RETURN
      END


