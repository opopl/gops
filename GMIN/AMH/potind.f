
C     --------------------- POTIND ----------------------

      SUBROUTINE POTIND(MAXSIZ,MAXCNT,NUMLNG,ILONG,
     *                  JBEGN,JEND,NMRES,I_TAB)

C     ---------------------------------------------------

C     POTIND SETS UP CONSTRAINT LIST FOR SPECIFIED
C            PARAMETERS (LSTRT,LFINS)

C     ARGUMENTS:

C        NMRES - ACTUAL NUMBER OF RESIDUES

C     ---------------------------------------------------

C     SET REQUIRED PARAMETERS

      USE AMHGLOBALS,  ONLY: IXN_FROM_SITE,NUMCONST_A,SEGLIST_A,MAXRES,
     *                I_IXN_QBIAS_B , I_IXN_QBIAS_A,
     *                NUMCONST_B, SEGLIST_B
      IMPLICIT NONE

C     ARGUMENT DECLARATIONS:

         INTEGER MAXSIZ,MAXCNT,NUMLNG(0:MAXSIZ),
     *           ILONG(MAXCNT,2),NMRES,
     *           JBEGN,JEND,I_TAB,
     *           ISIT1,ISIT2,I_TEMP   
  
C     INTERNAL VARIABLES:

         INTEGER IDD

C        --- DO LOOP INDICES ---

         INTEGER I_RES,I,J
 
C        --- IMPLIED DO LOOP INDICES ---

         INTEGER RES_A, RES_B,ARSEFLAP

C     --------------------- BEGIN -----------------------


C     SET UP INDICES FOR PAIRS OF CONSTRAINTS UTILIZED IN
C     POTENTIAL

      NUMLNG(0)=0
      IDD=0
      ARSEFLAP=0
      DO 500 RES_A=1,NMRES-JBEGN
         DO 901 RES_B=JBEGN,JEND-RES_A 
            ARSEFLAP=ARSEFLAP+1

C  ASSIGNMENT OF  ILONG

            ILONG(IDD+RES_B-JBEGN+1,1)=RES_A
            ILONG(IDD+RES_B-JBEGN+1,2)=RES_A + RES_B

            IXN_FROM_SITE(RES_A,RES_A+RES_B,I_TAB)
     *                             =IDD+RES_B-JBEGN+1
           

  901    CONTINUE 
         IDD=IDD + JEND - RES_A - JBEGN  + 1
         NUMLNG(RES_A)=IDD 

  500 CONTINUE

C          WRITE(6,*)'ARSEFLAP ',ARSEFLAP

      DO 502 I_RES=NMRES-JBEGN+1,NMRES
         NUMLNG(I_RES)=NUMLNG(NMRES-JBEGN)
  502 CONTINUE

C FINDS PROPER INTERACTION FROM LIST FOR BIASING OF SEGMENTS
C  CREAT TABLE ONCE

       IF (I_TAB .EQ. 1)THEN

       DO 30  I = 1, MAXRES
         DO 40 J = 1, MAXRES
                I_IXN_QBIAS_A(I,J)=0
40     CONTINUE
30     CONTINUE   

       DO 50 I = 1,  NUMCONST_A - 2
             ISIT1=SEGLIST_A(I)
        DO 100 J = I+2, NUMCONST_A
                ISIT2 = SEGLIST_A(J)
             DO I_TEMP=1,NUMLNG(NMRES)
                IF (ISIT1 .EQ. ILONG(I_TEMP,1) .AND.
     *                (ISIT2 ).EQ. ILONG(I_TEMP,2)) THEN
                        I_IXN_QBIAS_A(ISIT1,ISIT2) = I_TEMP
                ENDIF
             ENDDO ! I_TEMP=1,NUMLNG(NMRES)
100         CONTINUE
50          CONTINUE

       ENDIF ! IF (I_TAB .EQ. 1)THEN

C FINDS PROPER INTERACTION FROM LIST FOR BIASING OF SEGMENTS
C  CREAT TABLE ONCE

      IF (I_TAB .EQ. 1)THEN
       DO 301 I = 1, MAXRES
         DO 401 J = 1, MAXRES
             I_IXN_QBIAS_B(I,J)=0
401     CONTINUE
301     CONTINUE

       DO 501 I = 1,  NUMCONST_B - 2
             ISIT1=SEGLIST_B(I)
        DO 1001 J = I+2, NUMCONST_B
                ISIT2 = SEGLIST_B(J)
             DO I_TEMP=1,NUMLNG(NMRES)
                IF (ISIT1 .EQ. ILONG(I_TEMP,1) .AND.
     *                (ISIT2 ).EQ. ILONG(I_TEMP,2)) THEN
                    I_IXN_QBIAS_B(ISIT1,ISIT2) = I_TEMP
                ENDIF
             ENDDO ! I_TEMP=1,NUMLNG(NMRES)
1001     CONTINUE
501    CONTINUE
      ENDIF ! IF (I_TAB .EQ. 1)THEN

C     ---------------------- DONE -----------------------

      RETURN
      END
