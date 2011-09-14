C
C  Program gaussconv.f reads in a set of numbers and makes them
C  into a continuous distribution by convolution with Gaussians.
C
      PROGRAM GAUSSCONV

      IMPLICIT NONE

      INTEGER J1, J2, MMAX, N, NPTS, NEXTRA, NPOW, NXMAX, NXMIN, NVMAX, NVMIN
      PARAMETER(MMAX=10000,NPOW=10)
      DOUBLE PRECISION X(MMAX), DUMMY, XMAX, XMIN, DIST(MMAX), XVAL(MMAX), SIG, XNORM, PI, XSIG(NPOW), XMEAN(NPOW),
     1                  DUM(NPOW,MMAX), DUMMYN, SIGN, XNORMN, XVALN(MMAX), DISTN(MMAX), XMEANN(NPOW), XSIGN(NPOW),
     2                  DDUM(NPOW), DUMMYP, DUMMYPP, NXMEAN(NPOW), NXSIG(NPOW), DUMMYPN, DUMMYPPN, DUMN(NPOW,MMAX),
     3                  DDUMN(NPOW), NVMEAN(NPOW), NVSIG(NPOW), NX(MMAX), NV(MMAX)

      XMAX=-1.0D100
      XMIN=1.0D100
      NXMAX=-100000
      NXMIN= 100000
      NVMAX=-100000
      NVMIN= 100000
      XMEAN(1)=0.0D0
      XSIG(1)=0.0D0
      NXMEAN(1)=0.0D0
      NVMEAN(1)=0.0D0
      NXSIG(1)=0.0D0
      NVSIG(1)=0.0D0
      
      DO J1=1,MMAX
         READ(*,*,END=20) NX(J1), NV(J1), X(J1)
         XMEAN(1)=XMEAN(1)+ X(J1)
         XSIG(1)=XSIG(1)+ X(J1)**2
         IF (X(J1).GT.XMAX) XMAX=X(J1)
         IF (X(J1).LT.XMIN) XMIN=X(J1)
         NXMEAN(1)=NXMEAN(1)+ NX(J1)
         NVMEAN(1)=NVMEAN(1)+ NV(J1)
         NXSIG(1)=NXSIG(1)+ NX(J1)**2
         NVSIG(1)=NVSIG(1)+ NV(J1)**2
         IF (NX(J1).GT.NXMAX) NXMAX=NX(J1)
         IF (NX(J1).LT.NXMIN) NXMIN=NX(J1)
         IF (NV(J1).GT.NVMAX) NVMAX=NV(J1)
         IF (NV(J1).LT.NVMIN) NVMIN=NV(J1)
      ENDDO
20    N=J1-1
      XMEAN(1)=XMEAN(1)/N
      NXMEAN(1)=NXMEAN(1)/N
      NVMEAN(1)=NVMEAN(1)/N
      XSIG(1)=DSQRT((XSIG(1)-N*XMEAN(1)**2)/(N-1))
      NXSIG(1)=DSQRT((NXSIG(1)-N*NXMEAN(1)**2)/(N-1))
      NVSIG(1)=DSQRT((NVSIG(1)-N*NVMEAN(1)**2)/(N-1))
      WRITE(*,'(A,6F13.2)') 'mu/sig of t, Q, V  are ',
     1             XMEAN(1),XSIG(1),NXMEAN(1),NXSIG(1),NVMEAN(1),NVSIG(1)

      STOP

      NPTS=2000
      SIG=100.0D0
      XNORM=0.0D0
      SIGN=500.0D0
      XNORMN=0.0D0
      NEXTRA=20
      PI=3.141592654D0
      DO J1=1,NPTS+2*NEXTRA
         XVAL(J1)=XMIN+(XMAX-XMIN)*(J1-NEXTRA)/NPTS
         XVALN(J1)=NXMIN+(NXMAX-NXMIN)*(J1-NEXTRA)/NPTS
         DUMMYN=0.0D0
         DUMMY=0.0D0
         DO J2=1,N
            DUMMY=DUMMY+DEXP(-(XVAL(J1)-X(J2))**2/(2.0D0*SIG**2))/(DSQRT(2.0D0*PI)*SIG)
            DUMMYN=DUMMYN+DEXP(-(XVALN(J1)-NX(J2))**2/(2.0D0*SIGN**2))/(DSQRT(2.0D0*PI)*SIGN)
         ENDDO
         XNORM=XNORM+DUMMY
         XNORMN=XNORMN+DUMMYN
         DIST(J1)=DUMMY
         DISTN(J1)=DUMMYN
      ENDDO
      XNORM=XNORM*(XMAX-XMIN)/NPTS
      XNORMN=XNORMN*(NXMAX-NXMIN)/NPTS

      DUMMY=0.0D0
      DUMMYP=0.0D0
      DUMMYPP=0.0D0
      XMEAN(1)=0.0D0
      XSIG(1)=0.0D0
      DUMMYN=0.0D0
      DUMMYPN=0.0D0
      DUMMYPPN=0.0D0
      XMEANN(1)=0.0D0
      XSIGN(1)=0.0D0
      DO J1=2,NPOW
         DDUM(NPOW)=0.0D0
         XMEAN(J1)=0.0D0
         XSIG(J1)=0.0D0
         DDUMN(NPOW)=0.0D0
         XMEANN(J1)=0.0D0
         XSIGN(J1)=0.0D0
      ENDDO
      DO J1=1,NPTS+2*NEXTRA
         DIST(J1)=DIST(J1)/XNORM
         DUMMY=DUMMY+DIST(J1)*(XMAX-XMIN)/NPTS
         XMEAN(1)=XMEAN(1)+XVAL(J1)*DIST(J1)*(XMAX-XMIN)/NPTS
         DISTN(J1)=DISTN(J1)/XNORMN
         DUMMYN=DUMMYN+DISTN(J1)*(NXMAX-NXMIN)/NPTS
         XMEANN(1)=XMEANN(1)+XVALN(J1)*DISTN(J1)*(NXMAX-NXMIN)/NPTS
         DO J2=2,NPOW
            DUM(J2,J1)=  (1-DUMMYP)**(J2-1)*DIST(J1)*J2
            DUMN(J2,J1)=  (1-DUMMYPN)**(J2-1)*DISTN(J1)*J2
         ENDDO
         WRITE(*,'(4F20.10)') XVAL(J1),DIST(J1),XVALN(J1),DISTN(J1)
         DUMMYPP=DUMMYP
         DUMMYP=DUMMY
         DUMMYPPN=DUMMYPN
         DUMMYPN=DUMMYN
         DO J2=2,NPOW
            DDUMN(J2)=DDUMN(J2)+DUMN(J2,J1)*(NXMAX-NXMIN)/NPTS
         ENDDO
      ENDDO
C     PRINT*,'DUMMY=',DUMMY
C     PRINT*,'DDUM=',DDUM
      DO J2=2,10
         DO J1=1,NPTS+2*NEXTRA
C           WRITE(*,'(2F20.10)') J2*XVAL(J1),DUM(J2,J1)/J2
            XMEAN(J2)=XMEAN(J2)+J2*XVAL(J1)*DUM(J2,J1)*(XMAX-XMIN)/NPTS
            XMEANN(J2)=XMEANN(J2)+J2*XVALN(J1)*DUMN(J2,J1)*(NXMAX-NXMIN)/NPTS
         ENDDO
      ENDDO
      DO J1=1,NPTS+2*NEXTRA
         XSIG(1)=XSIG(1)+(XVAL(J1)-XMEAN(1))**2*DIST(J1)*(XMAX-XMIN)/NPTS
         XSIGN(1)=XSIGN(1)+(XVALN(J1)-XMEANN(1))**2*DISTN(J1)*(NXMAX-NXMIN)/NPTS
         DO J2=2,NPOW
            XSIG(J2)=XSIG(J2)+(J2*XVAL(J1)-XMEAN(J2))**2*DUM(J2,J1)*(XMAX-XMIN)/NPTS
            XSIGN(J2)=XSIGN(J2)+(J2*XVALN(J1)-XMEANN(J2))**2*DUMN(J2,J1)*(NXMAX-NXMIN)/NPTS
         ENDDO
      ENDDO
      DDUM(1)=DUMMY
      WRITE(*,'(A,10F10.3)') '# Normalization=',DUMMY,(DDUM(J1),J1=2,NPOW)
      WRITE(*,'(A,10F10.3)') '# Mean         =',(XMEAN(J1)/DDUM(J1),J1=1,NPOW)
      WRITE(*,'(A,10F10.3)') '# Sigma        =',(DSQRT(XSIG(J1))/DDUM(J1),J1=1,NPOW)

      STOP
      END
