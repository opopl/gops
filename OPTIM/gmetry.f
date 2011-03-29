C   OPTIM: A PROGRAM FOR OPTIMIZING GEOMETRIES AND CALCULATING REACTION PATHWAYS
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF OPTIM.
C
C   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
C
C     ROUTINE TO CALCULATE CARTESIAN COORDINATES FROM INTERNAL
C     COORDINATE REPRESENTATION.  SOME OF THIS HAS BEEN LIFTED FROM
C     PRDDO, ALTHOUGH SOME IMPROVEMENTS HAVE BEEN MADE.
C     CONNECTIVITY OF FIRST THREE MUST BE 1-2-3 IN
C     INTERNAL COORDINATE REP.
C
      SUBROUTINE GMETRY(ITER,VEC,Q)
      USE KEY
      USE COMMONS
      IMPLICIT NONE
      INTEGER ITER, NUMAT, NDIS, I, J1, J2, J, NCOUNT
      DOUBLE PRECISION VEC(3*NATOMS), Q(3*NATOMS), SCALE, TX, TY, TZ
C
      NUMAT=0
      NDIS=1
      SCALE=1.0D0
      IF (IPRNT .GE. 20) THEN
         WRITE (*, *) 'CARTESIAN COORDINATES BEFORE SCALING.'
         WRITE (*, '(I3,3F12.6)') (I,Q(3*I-2),Q(3*I-1),Q(3*I),
     1      I=1,NATOMS)
      ENDIF
C     IF (NDIS.NE.0) SCALE=0.5291772D0
      IF (NDIS.NE.0) SCALE=1.0D0
      DO I=1,NATOMS
         Q(3*I-2)=Q(3*I-2)/SCALE
         Q(3*I-1)=Q(3*I-1)/SCALE
         Q(3*I)=Q(3*I)/SCALE
      ENDDO
      IF (IPRNT .GE. 20) THEN
         WRITE (*, *) 'CARTESIAN COORDINATES AFTER SCALING.'
         WRITE (*, '(I3,3F12.6)') (I,Q(3*I-2),Q(3*I-1),Q(3*I),
     1      I=1,NATOMS)
      ENDIF
C
C PERMUTE THE COORDINATES IF NECESSARY - ONLY FOR FH
C
      IF (ZSYM(NATOMS)(1:2).EQ.'FH') THEN
      NCOUNT=0
130   CONTINUE
      DO 150 J1=1,NATOMS-1
         TX=Q(3*(J1-1)+1)-Q(3*J1+1)
         TY=Q(3*(J1-1)+2)-Q(3*J1+2)
         TZ=Q(3*(J1-1)+3)-Q(3*J1+3)
         IF (DABS(TZ*TZ+TY*TY).LT.1.0D-10) THEN
            NCOUNT=NCOUNT+1
            IF (NCOUNT.GT.2) THEN
               PRINT*,'*** WARNING ***'
               PRINT*,'*** CANNOT PERMUTE THE COORDINATES IN GMETRY'
               GOTO 160
C              STOP
            ENDIF
            DO 140 J2=1,NATOMS
               TX=Q(3*(J2-1)+1)
               TY=Q(3*(J2-1)+2)
               TZ=Q(3*(J2-1)+3)
               Q(3*(J2-1)+1)=TY
               Q(3*(J2-1)+2)=TZ
               Q(3*(J2-1)+3)=TX
140         CONTINUE
            GOTO 130
         ENDIF
150   CONTINUE
      ENDIF

160   IF (NATOMS.EQ.NUMAT) RETURN
      J=0
      DO 170 I=1,NATOMS
         IF (NR(I) .NE. 99) THEN
            J=J+1
C
C  TURN OFF THE AXIS PERMUTATIONS -
C  THIS IS THE ONLY WAY I CAN FIND OF PRESERVING THE
C  ORIGINAL CARTESIAN COORDINATES.
C
            Q(3*J-2)=Q(3*I-2)
            Q(3*J-1)=Q(3*I-1)
            Q(3*J)=Q(3*I)
         ENDIF
170   CONTINUE
      IF (IPRNT .GE. 20) THEN
         WRITE (*,*) 'CARTESIAN COORDINATES FROM GMETRY.'
C
         WRITE (*,'(I3,3F12.6)') (I,Q(3*I-2),Q(3*I-1),Q(3*I),I=1,NATOMS)
      ENDIF
C
C CHECK THE THETA ANGLES IF WE ARE DOING WATER.
C
      IF ((ZSYM(NATOMS)(1:1).EQ.'W').AND.(.NOT.ANGLEAXIS)) CALL TCHECK(VEC,Q)
      RETURN
      END
