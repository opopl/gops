
      SUBROUTINE WALESAMH_INTERFACE(COORD_MCP,GRAD_FOR_WALES,E_FOR_WALES)

!     CALCULATES ENERGIES AND FORCES FOR GMIN VIA A SCLTAB LOOKUP

     USE AMHGLOBALS, ONLY: MAXTAB,MAXR,PRCORD,ZRCORD,AVEP,VPOTNT, &
       FORSE,AMHMAXSIZ,NMRES,PEXCLD,NUMLNG,RINCINV, RINCSQ,ILONG,CRDIXN, &
       IRES,EQDIST,HBOND,OXEXCLDV,AVEP,MAXPRO,MAXCRD,X_MCP, WEIGHT_P_AP

      IMPLICIT NONE

!     INTERNAL VARIABLES:

      LOGICAL TEMPAV,SCL_CALL
      INTEGER JSTRT,JFINS,GLY_COUNT,III

      DOUBLE PRECISION GRAD_FOR_WALES(NMRES*3*3),E_FOR_WALES,COORD_MCP(NMRES*3*3)
      DOUBLE PRECISION F_CORD(AMHMAXSIZ,3,MAXCRD),E_TEMP_MCP,TRGENG(MAXTAB,3)
      DOUBLE PRECISION E_P_AP,E_PE_NO_BIAS

!     REQUIRED SUBROUTINES:
      EXTERNAL FORCE, HARM_SPRING

!     --------------------- BEGIN -----------------------
!  FIND POTENTIAL ENERGY FOR EACH INTERACTION TYPE

      TEMPAV=.TRUE.
      GLY_COUNT = 0

      DO III = 1,NMRES
       IF (IRES(III).EQ.8) THEN
        PRCORD(III, 1, 1, 1) = (COORD_MCP(9*(III-1)+1- GLY_COUNT*3)) !  CA X
        PRCORD(III, 2, 1, 1) = (COORD_MCP(9*(III-1)+2- GLY_COUNT*3)) !  CA Y
        PRCORD(III, 3, 1, 1) = (COORD_MCP(9*(III-1)+3- GLY_COUNT*3)) !  CA Z
!    SWAP  CA FOR CB
        PRCORD(III, 1, 1, 2) = (COORD_MCP(9*(III-1)+1- GLY_COUNT*3)) !  CB X
        PRCORD(III, 2, 1, 2) = (COORD_MCP(9*(III-1)+2- GLY_COUNT*3)) !  CB Y
        PRCORD(III, 3, 1, 2) = (COORD_MCP(9*(III-1)+3- GLY_COUNT*3)) !  CB Z
        PRCORD(III, 1, 1, 3) = (COORD_MCP(9*(III-1)+4- GLY_COUNT*3)) !  O X
        PRCORD(III, 2, 1, 3) = (COORD_MCP(9*(III-1)+5- GLY_COUNT*3)) !  O Y
        PRCORD(III, 3, 1, 3) = (COORD_MCP(9*(III-1)+6- GLY_COUNT*3)) !  O Z
        GLY_COUNT = GLY_COUNT +1
       ELSE
        PRCORD(III, 1, 1, 1) = (COORD_MCP(9*(III-1)+1- GLY_COUNT*3)) !  CA X
        PRCORD(III, 2, 1, 1) = (COORD_MCP(9*(III-1)+2- GLY_COUNT*3)) !  CA Y
        PRCORD(III, 3, 1, 1) = (COORD_MCP(9*(III-1)+3- GLY_COUNT*3)) !  CA Z
        PRCORD(III, 1, 1, 2) = (COORD_MCP(9*(III-1)+4- GLY_COUNT*3)) !  CB X
        PRCORD(III, 2, 1, 2) = (COORD_MCP(9*(III-1)+5- GLY_COUNT*3)) !  CB Y
        PRCORD(III, 3, 1, 2) = (COORD_MCP(9*(III-1)+6- GLY_COUNT*3)) !  CB Z
        PRCORD(III, 1, 1, 3) = (COORD_MCP(9*(III-1)+7- GLY_COUNT*3)) !  O X
        PRCORD(III, 2, 1, 3) = (COORD_MCP(9*(III-1)+8- GLY_COUNT*3)) !  O Y
        PRCORD(III, 3, 1, 3) = (COORD_MCP(9*(III-1)+9- GLY_COUNT*3)) !  O Z
       ENDIF
      ENDDO       

      SCL_CALL=.TRUE.

      JSTRT=1
      JFINS=NMRES

        E_TEMP_MCP=0.0D0
        E_P_AP=0.0D0           
        E_FOR_WALES=0.0D0 
        E_PE_NO_BIAS=0.0D0
        ZRCORD(:,:,1,:)=0.D0

       CALL HARM_SPRING(PRCORD,JSTRT,JFINS,MAXPRO,MAXCRD,IRES,F_CORD,E_TEMP_MCP)

       CALL FORCE(1,PRCORD,ZRCORD,AVEP,TEMPAV,MAXR,VPOTNT,FORSE,TRGENG,PEXCLD, &
              NUMLNG,NMRES,RINCINV,RINCSQ,ILONG,CRDIXN,IRES,EQDIST,HBOND,OXEXCLDV,SCL_CALL)

        ZRCORD(:,:,1,:)=ZRCORD(:,:,1,:)+F_CORD

         E_FOR_WALES =                  & !
                        AVEP(1,1,1)  +   & ! HYDRN POTENTIAL
                        AVEP(1,1,2)  +   & ! RAMA POTENTIAL
                        AVEP(1,1,3)  +   & ! PEPTIDE BOND
                        AVEP(1,1,4)  +   & ! CHIRALITY
                        AVEP(1,1,5)  +   & ! TOTAL AMH
!                       AVEP(1,1,6)  +   & ! ?????  EMPTY
!                       AVEP(1,1,7)  +   & ! AMHSR
!                       AVEP(1,1,8)  +   & ! AMHLR
                        AVEP(1,1,9)  +   & ! CARBON EX VOL
!                       AVEP(1,1,10) +   & ! Q_BIAS_A
                        AVEP(1,1,11) +   & ! OXY EX VOL
!                       AVEP(1,1,12) +   & ! OBIAS_RG
!                       AVEP(1,1,13) +   & ! AMHMR
!                       AVEP(1,1,14) +   & ! OHDRGN_S   AVEP(I_PRO,1,21:24)
!                       AVEP(1,1,15) +   & ! OHDRGN_M   AVEP(I_PRO,1,25:28) 
!                       AVEP(1,1,16) +   & ! OHDRGN_L   AVEP(I_PRO,1,29:32)
                        AVEP(1,1,17) +   & ! NON_ADD  
!                       AVEP(1,1,18) +   & ! Q_BIAS_B
!                       AVEP(1,1,19) +   & ! ?????  EMPTY
!                       AVEP(1,1,20) +   & ! ?????  EMPTY
!                       AVEP(1,1,40) +   & ! E_P_AP_1
!                       AVEP(1,1,41) +   & ! E_P_AP_2
!                       AVEP(1,1,42) +   & ! E_P_AP_3
!                       E_P_AP       +   & ! E_P_AP
!                       AVEP(1,1,43) +   & ! REPLICA ENERGY TERM
                        E_TEMP_MCP

!     I_PRO                   INDEX OVER PROTEINS
!     I_STEP                  INDEX OVER TIME STEPS PER TEMPERATURE
!     I_TEMP                  INDEX OVER TEMPERATURES
!     ITCNT                   INCREMENT T (TEMPERATURE) COUNTER

        GLY_COUNT = 0
      DO III = 1,NMRES
       IF (IRES(III).EQ.8) THEN
         COORD_MCP(9*(III-1)+1- GLY_COUNT*3) = (PRCORD(III, 1, 1, 1)) !  CA X
         COORD_MCP(9*(III-1)+2- GLY_COUNT*3) = (PRCORD(III, 2, 1, 1)) !  CA Y
         COORD_MCP(9*(III-1)+3- GLY_COUNT*3) = (PRCORD(III, 3, 1, 1)) !  CA Z
         COORD_MCP(9*(III-1)+4- GLY_COUNT*3) = (PRCORD(III, 1, 1, 3)) !  O X
         COORD_MCP(9*(III-1)+5- GLY_COUNT*3) = (PRCORD(III, 2, 1, 3)) !  O Y
         COORD_MCP(9*(III-1)+6- GLY_COUNT*3) = (PRCORD(III, 3, 1, 3)) !  O Z
         GLY_COUNT = GLY_COUNT +1 
       ELSE
         COORD_MCP(9*(III-1)+1 - GLY_COUNT*3) = (PRCORD(III, 1, 1, 1)) !  CA X
         COORD_MCP(9*(III-1)+2 - GLY_COUNT*3) = (PRCORD(III, 2, 1, 1)) !  CA Y
         COORD_MCP(9*(III-1)+3 - GLY_COUNT*3) = (PRCORD(III, 3, 1, 1)) !  CA Z
         COORD_MCP(9*(III-1)+4 - GLY_COUNT*3) = (PRCORD(III, 1, 1, 2)) !  CB X
         COORD_MCP(9*(III-1)+5 - GLY_COUNT*3) = (PRCORD(III, 2, 1, 2)) !  CB Y
         COORD_MCP(9*(III-1)+6 - GLY_COUNT*3) = (PRCORD(III, 3, 1, 2)) !  CB Z
         COORD_MCP(9*(III-1)+7 - GLY_COUNT*3) = (PRCORD(III, 1, 1, 3)) !  O X
         COORD_MCP(9*(III-1)+8 - GLY_COUNT*3) = (PRCORD(III, 2, 1, 3)) !  O Y
         COORD_MCP(9*(III-1)+9 - GLY_COUNT*3) = (PRCORD(III, 3, 1, 3)) !  O Z
       ENDIF
      ENDDO

!       WRITE(6,'(5F10.5,2X)')E_PE_BACKBONE,AVEP(1,1,2),AVEP(1,1,5),E_P_AP,E_PE_NO_BIAS
!       WRITE(6,'(4F10.5,2X)')AVEP(1,1,1),AVEP(1,1,14),AVEP(1,1,15),AVEP(1,1,16)

!  AVEP(1,1,1,ITCNT)  = HYDROGEN BOND POTENTIAL
!  AVEP(1,1,2,ITCNT)  = RAMA POTENTIAL
!  AVEP(1,1,3,ITCNT)  = OXYGEN EXCLUDED
!  AVEP(1,1,4,ITCNT)  = CHIRALITY POTENTIAL
!  AVEP(1,1,5,ITCNT)  = TOTAL AMH SHORT + MEDIUM + LONG
!  AVEP(1,1,6,ITCNT)  = ------
!  AVEP(1,1,7,ITCNT)  = AMH SHORT RANGE
!  AVEP(1,1,8,ITCNT)  = AMH LONG RANGE
!  AVEP(1,1,9,ITCNT)  = POTENTIAL DUE TO CARBON (A/B) EXCLUDED VOLUME
!  AVEP(1,1,10,ITCNT) = Q BIASING POTENTAL
!  AVEP(1,1,11,ITCNT) = OXY POTENTIAL
!  AVEP(1,1,12,ITCNT) = RADIUS OF GYRATION BIAS
!  AVEP(1,1,13,ITCNT) = AMH MEDIUM RANGE
!  AVEP(1,1,14,ITCNT) = HYDROGEN BONDS SHORT RANGE
!  AVEP(1,1,15,ITCNT) = HYDROGEN BONDS MEDIUM RANGE
!  AVEP(1,1,16,ITCNT) = HYDROGEN BONDS LONG RANGE
!  AVEP(1,1,17,ITCNT) = NONADDATIVE
!  AVEP(1,1,18,ITCNT) = Q BIASING POTENTAL SEG

!        WRITE(6,*)'E SHAKE SPRINGS       ', E_TEMP_MCP
!        WRITE(6,*)'E 1  OHDRGN           ',AVEP(1,1,1)
!        WRITE(6,*)'E 2  ORAMA            ',AVEP(1,1,2)
!        WRITE(6,*)'E 3  OOXY             ',AVEP(1,1,3)
!        WRITE(6,*)'E 4  OCHIRAL          ',AVEP(1,1,4)
!        WRITE(6,*)'E 5  OAMH             ',AVEP(1,1,5)
!        WRITE(6,*)'E 6  ----             ',AVEP(1,1,6)
!        WRITE(6,*)'E 7  OAMHSR           ',AVEP(1,1,7)
!        WRITE(6,*)'E 8  OAMHLR           ',AVEP(1,1,8)
!        WRITE(6,*)'E 9  OCCEV            ',AVEP(1,1,9)
!        WRITE(6,*)'E 10 QBIAS_A          ',AVEP(1,1,10)
!        WRITE(6,*)'E 11 OOEV             ',AVEP(1,1,11)
!        WRITE(6,*)'E 12 RG_BIAS          ',AVEP(1,1,12)
!        WRITE(6,*)'E 13 AMHMR            ',AVEP(1,1,13)
!        WRITE(6,*)'E 14 HDRGNS           ',AVEP(1,1,14)
!        WRITE(6,*)'E 15 HDRGNM           ',AVEP(1,1,15)
!        WRITE(6,*)'E 16 HDRGNL           ',AVEP(1,1,16)
!        WRITE(6,*)'E 17 NON ADD CONTACT  ',AVEP(1,1,17)
!        WRITE(6,*)'E 18 QBIAS_B          ',AVEP(1,1,18)
!        WRITE(6,*)'E E_P_AP              ',E_P_AP 

        GLY_COUNT = 0

      DO III=1,NMRES
       IF (IRES(III).EQ.8) THEN
         GRAD_FOR_WALES(9*(III-1)+1 - GLY_COUNT*3) = -(ZRCORD(III,1,1,1)) ! CAX
         GRAD_FOR_WALES(9*(III-1)+2 - GLY_COUNT*3) = -(ZRCORD(III,2,1,1)) ! CAY 
         GRAD_FOR_WALES(9*(III-1)+3 - GLY_COUNT*3) = -(ZRCORD(III,3,1,1)) ! CAZ 
         GRAD_FOR_WALES(9*(III-1)+4 - GLY_COUNT*3) = -(ZRCORD(III,1,1,3)) ! OX
         GRAD_FOR_WALES(9*(III-1)+5 - GLY_COUNT*3) = -(ZRCORD(III,2,1,3)) ! OY 
         GRAD_FOR_WALES(9*(III-1)+6 - GLY_COUNT*3) = -(ZRCORD(III,3,1,3)) ! OZ 
         GLY_COUNT = GLY_COUNT +1
       ELSE  
         GRAD_FOR_WALES(9*(III-1)+1 - GLY_COUNT*3) = -(ZRCORD(III,1,1,1)) ! CAX
         GRAD_FOR_WALES(9*(III-1)+2 - GLY_COUNT*3) = -(ZRCORD(III,2,1,1)) ! CAY 
         GRAD_FOR_WALES(9*(III-1)+3 - GLY_COUNT*3) = -(ZRCORD(III,3,1,1)) ! CAZ 
         GRAD_FOR_WALES(9*(III-1)+4 - GLY_COUNT*3) = -(ZRCORD(III,1,1,2)) ! CBX 
         GRAD_FOR_WALES(9*(III-1)+5 - GLY_COUNT*3) = -(ZRCORD(III,2,1,2)) ! CBY 
         GRAD_FOR_WALES(9*(III-1)+6 - GLY_COUNT*3) = -(ZRCORD(III,3,1,2)) ! CBZ 
         GRAD_FOR_WALES(9*(III-1)+7 - GLY_COUNT*3) = -(ZRCORD(III,1,1,3)) ! OX
         GRAD_FOR_WALES(9*(III-1)+8 - GLY_COUNT*3) = -(ZRCORD(III,2,1,3)) ! OY 
         GRAD_FOR_WALES(9*(III-1)+9 - GLY_COUNT*3) = -(ZRCORD(III,3,1,3)) ! OZ 
       ENDIF
      ENDDO

      RETURN
      END
