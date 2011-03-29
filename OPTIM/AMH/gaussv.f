
C     --------------------- GAUSSV ----------------------

      SUBROUTINE GAUSSV(MAXS,DIST,DELTR,GAUSSP,RGRID)

C     ---------------------------------------------------

C     GAUSSV GENERATES A GUASSIAN CURVE CENTERED ABOUT 
C            DIST WITH A 'STANDARD DEVIATION' OF DELTR

C     ARGUMENTS:

C        MAXS  - NUMBER OF POINTS THE GRID IS TO BE
C                EVALUATED AT (I)
C        DIST  - CENTER OR AVERAGE OF GAUSSIAN; 
C                IT IS THE DISTANCE BETWEEN RESIDUES 
C                I AND J (I)
C        DELTR - IS THE WIDTH OR STANDARD DEVIATION OF
C                THE WELL (I)
C        GAUSSP- GAUSSIAN AS A FUNCTION OF THE R-GRID
C                IN RGRID (O)
C        RGRID - GRID OF R POINTS FOR WHICH THE GAUSSIAN
C                IS TO BE COMPUTED (I)        

C     ---------------------------------------------------

      IMPLICIT NONE

C     ARGUMENT DECLARATIONS:

         INTEGER MAXS

         DOUBLE PRECISION DIST,DELTR,GAUSSP(MAXS),
     *        RGRID(MAXS)


C     INTERNAL VARIABLES:
 
         INTEGER I500,I501,I502,I503

C     --------------------- BEGIN -----------------------

C     --- DIAGNOSTICS ---

C     ECHO INPUT
 
C      WRITE(OARCHV,100)MAXS,DIST,DELTR
C  100 FORMAT('GAUSSV:MAXS ',I3,' DIST AND DELTR',
C    *        2(1X,1PE10.3))

C     --- END DIAGNOSTICS ---


C     FIND (R - R(MEM))**2

      DO 500 I500=1,MAXS
         GAUSSP(I500)=-( RGRID(I500) - DIST )**2
  500 CONTINUE

C     DIVIDE BY 'VARIANCE' OF WELL-WIDTH
            
      DO 501 I501=1,MAXS
         GAUSSP(I501)=GAUSSP(I501)*DELTR
  501 CONTINUE

C     MAKE SURE UNDERFLOW DOESN'T OCCUR

      DO 502 I502=1,MAXS
         GAUSSP(I502)=MAX( GAUSSP(I502),-60.0D0 )
  502 CONTINUE

C     COMPUTE GAUSSIAN

      DO 503 I503=1,MAXS
         GAUSSP(I503)=EXP( GAUSSP(I503) )
  503 CONTINUE




C     ---------------------- DONE -----------------------
 
      RETURN
      END
