      SUBROUTINE WALESAMH_INITIAL()

      USE AMHGLOBALS,  ONLY:SO, NMTEMP,ITGRD,TEMGRD,TEMTUR,ICTEMP,CTEMP,
     *  ISCLTAB,NMRES,OARCHV,NMSTEP,NUMPRO,IDIGNS,MAXPRO,MAXCRD,PRCORD,IRES,OCONV,
     *  OMOVI,OMOVISEG,QUENCH,NQUENCH,QUENCH_CRD

      IMPLICIT NONE

C     SUBROUTINES REQUIRED BY MAIN PROGRAM

       EXTERNAL GENTAB,INITIL,INTSTR,SCLTAB,ZERO_AMH

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     INTERNAL VARIABLES:

         INTEGER JSTRT,JFINS,I_QUENCH,LEN,ISHKIT,NMDIFV

         CHARACTER*10 SAVE_NAME

      CALL ZERO_AMH

C     --------------------- BEGIN -----------------------


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     OPEN REQUIRED FILES, READ INPUT PARAMETER FILE, AND GENERATE HEADER FILE
C      CALL READ_INPUT_ALT() ! CALLED BEFORE INITIL : JOHAN
      CALL INITIL
C      CALL READ_ALTGAMMA()  ! CALLED AFTER  INITIL : JOHAN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SET UP TEMPERATURE-ANNEALING SCHEDULE
C      CALL ANNSCH(NMTEMP,ITGRD,TEMGRD,TEMTUR,ICTEMP,CTEMP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     GENERATE REQUISITE FORCE/POTENTIAL TABLES
      CALL GENTAB
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     SCALE TABLES
      IF (ISCLTAB) WRITE(SO,*) 'IN SCLTAB'
      IF (ISCLTAB) CALL SCLTAB
      IF (ISCLTAB) WRITE(SO,*) 'OUT SCLTAB'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     GENERATE INITIAL STRUCTURES
      QUENCH_CRD=0.0D0
    
C      WRITE(6,*) 'IN INTSTR'
      CALL INTSTR
C      WRITE(6,*) 'OUT INTSTR'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        
        IF (.NOT. QUENCH) NQUENCH=1
        DO I_QUENCH = 1,NQUENCH
        PRCORD=QUENCH_CRD(:,:,:,:,I_QUENCH)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SET INDICIES FOR THE FIRST AND LAST RESIDUES
C     WHICH ARE NOT FIXED IN CRYSTAL CONFORMATION
      JSTRT=1
      JFINS=NMRES

C     SET SUBSEGMENT LENGTH
      LEN=JFINS - JSTRT + 1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C++++++++++++++++++++++++++++++JOHAN
!       CALL READ_INPUT_ALT()
!       CALL READ_ALTGAMMA() ! MUST BE CALLED AFTER INITIL

C------------------------------JOHAN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     --- DIAGNOSTICS ---
C      WRITE(OARCHV,121)JSTRT,JFINS,LEN
C  121 FORMAT(/' START ',I3,' END ',I3,' LENGTH ',I3)
C      WRITE(OARCHV,122)JSTRT,JFINS,NMSTEP,NUMPRO
C  122 FORMAT('JSTRT ',I3,' JFINS ',I3,' MUTATIONS/T ',
C     *       I3,' NUMPRO ',I3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     GENERATE INITIAL ENSEMBLE OF PROTEINS
C     FIND CONFIGURATION WHICH SATISFIES THE CONSTRAINTS

      IDIGNS=.FALSE.

        ENDDO ! I_QUENCH

      END
