      SUBROUTINE EV_SET_UP

      USE AMHGLOBALS,  ONLY:SO, EXVMIN_GAMMA,NMRES,IRES,NUMLNG,EXVMINS,EXVMIN,&
            CCEV_DIST,ILONG,CRDIXN,IEXCLD_BETA,IEXCLD_GAMMA,C_OF_M_DIST,&
            EXVMINS_BETA,O_EXVMINS,O_EXVMIN,OOEV_DIST,IEVGAMMA,IEVBETA

      IMPLICIT NONE

      DOUBLE PRECISION HRDRAD_GAMMA(1:20),HRDRAD_BETA(1:20)
       DOUBLE PRECISION, PARAMETER, DIMENSION(20):: C_OF_M_IN_ANGSTROMS=(/1.52D0,3.56D0,2.09D0,2.09D0,1.88D0,&
                                                             2.63D0,2.69D0,0.00D0,2.39D0,2.19D0,&
                                                             2.31D0,2.99D0,2.55D0,2.52D0,1.42D0,&
                                                             1.78D0,1.89D0,2.73D0,2.67D0,1.91D0/)

      INTEGER I,J,I_TAB,I_IXN,I_ATOM,J_ATOM,OPEN_STATUS,I_DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CENTRE OF MASS DISTANCES FOR `GAMMA' EXCLUDED VOLUME

      DO I=1,NMRES   !IN UNITS OF CA-CB BOND LENGTH
        C_OF_M_DIST(I)=C_OF_M_IN_ANGSTROMS(IRES(I))/1.52D0
      ENDDO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GAMMA EXCLUDED VOLUME

      EXVMIN_GAMMA=0.0D0
      
      IF (IEXCLD_GAMMA) THEN
         OPEN(UNIT=IEVGAMMA,FILE='PARAMS/R_EV_GAMMA.DAT',&
           STATUS='OLD',ACTION='READ',IOSTAT=OPEN_STATUS)
         IF (OPEN_STATUS.NE.0) THEN
             WRITE(SO,*) 'ERROR OPENING GAMMA EV FILE'
             STOP
         ENDIF
         DO I=1,20
           READ(IEVGAMMA,*) I_DUMMY,HRDRAD_GAMMA(I)
         ENDDO
      ENDIF

      DO I=1,NMRES-13
      DO J=I+13,NMRES
        EXVMIN_GAMMA(I,J)=HRDRAD_GAMMA(IRES(I))+HRDRAD_GAMMA(IRES(J))
      ENDDO
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (IEXCLD_BETA) THEN !SET *BETA* VALUES. ONLY USED IF IEXCLD_BETA=.TRUE.
         OPEN(UNIT=IEVBETA,FILE='PARAMS/R_EV_BETA.DAT',&
           STATUS='OLD',ACTION='READ',IOSTAT=OPEN_STATUS)
         IF (OPEN_STATUS.NE.0) THEN
             WRITE(SO,*) 'ERROR OPENING BETA EV FILE'
             STOP
         ENDIF
         DO I=1,20
           READ(IEVBETA,*) I_DUMMY,HRDRAD_BETA(I)
         ENDDO
      ENDIF

         DO I_TAB=1,3             !LEAVE BETA-BETA AS SPECIAL CASE
         DO I_IXN=1,NUMLNG(NMRES,I_TAB)
            I=ILONG(I_IXN,1,I_TAB)
            J=ILONG(I_IXN,2,I_TAB)

            I_ATOM=CRDIXN(I_TAB,1)
            J_ATOM=CRDIXN(I_TAB,2)

            IF( (I_ATOM.EQ.2).AND.(IRES(I).EQ.8) ) CYCLE
            IF( (J_ATOM.EQ.2).AND.(IRES(J).EQ.8) ) CYCLE

!              SET EQUILIBRIUM DISTANCE

            IF (IABS(J-I).LT.5) THEN
              CCEV_DIST(I_IXN,I_TAB)=EXVMINS 
            ELSE
              CCEV_DIST(I_IXN,I_TAB)=EXVMIN
            ENDIF

          ENDDO
          ENDDO

          I_TAB=4
          DO I_IXN=1,NUMLNG(NMRES,I_TAB)
            I=ILONG(I_IXN,1,I_TAB)
            J=ILONG(I_IXN,2,I_TAB)
            I_ATOM=CRDIXN(I_TAB,1)
            J_ATOM=CRDIXN(I_TAB,2)

            IF( (I_ATOM.EQ.2).AND.(IRES(I).EQ.8) ) CYCLE
            IF( (J_ATOM.EQ.2).AND.(IRES(J).EQ.8) ) CYCLE

!              SET EQUILIBRIUM DISTANCE

            IF (IABS(J-I).LT.5) THEN
              IF (IABS(J-I).LT.1) THEN
                WRITE(SO,*) 'PROBLEM SETTING UP EV',IABS(J-I),I_IXN
                STOP
              ENDIF
              CCEV_DIST(I_IXN,I_TAB)=EXVMINS_BETA(IABS(J-I)) 
            ELSEIF (IABS(J-I).LT.13) THEN
              CCEV_DIST(I_IXN,I_TAB)=EXVMIN
            ELSEIF(.NOT.IEXCLD_BETA) THEN
              CCEV_DIST(I_IXN,I_TAB)=EXVMIN
            ELSE 
              CCEV_DIST(I_IXN,I_TAB)=HRDRAD_BETA(IRES(I))+HRDRAD_BETA(IRES(J))
            ENDIF
          ENDDO

! DO O-O EXCLUDED VOLUME

          I_TAB=2
          DO I_IXN=1,NUMLNG(NMRES,I_TAB)

            I=ILONG(I_IXN,1,I_TAB)
            J=ILONG(I_IXN,2,I_TAB)

            IF (IABS(J-I).LT.5) THEN
              OOEV_DIST(I_IXN)=O_EXVMINS(IABS(J-I))
            ELSE
              OOEV_DIST(I_IXN)=O_EXVMIN
            ENDIF
          ENDDO

          WRITE(SO,*) '*****************************'
          WRITE(SO,*) 'O_EXVMIN=',O_EXVMIN
          WRITE(SO,*) 'O_EXVMINS=',O_EXVMINS
          WRITE(SO,*) '*****************************'
          WRITE(SO,*) 'EXVMIN=',EXVMIN
          WRITE(SO,*) 'EXVMINS=',EXVMINS
          WRITE(SO,*) 'EXVMINS_BETA=',EXVMINS_BETA
          WRITE(SO,*) '*****************************'

          RETURN
          END
      
