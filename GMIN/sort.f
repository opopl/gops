C   GMIN: A PROGRAM FOR FINDING GLOBAL MINIMA
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF GMIN.
C
C   GMIN IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   GMIN IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
      SUBROUTINE GSORT(N,ARR,ARR2,NP,NPAR,NTAB)
      INTEGER N,M,NSTACK,NP,NPAR,NTAB,A2
      DOUBLE PRECISION ARR(NTAB,NPAR), ARR2(NTAB,NPAR)
      PARAMETER (M=7,NSTACK=50)
      INTEGER I,IR,J,JSTACK,K,L,ISTACK(NSTACK)
      DOUBLE PRECISION A,TEMP,TEMP2
      JSTACK=0
      L=1
      IR=N
1     IF(IR-L.LT.M)THEN
        DO 12 J=L+1,IR
          A=ARR(J,NP)
C******************************
          A2=ARR2(J,NP)
C******************************
          DO 11 I=J-1,1,-1
            IF(ARR(I,NP).LE.A)GOTO 2
            ARR(I+1,NP)=ARR(I,NP)
C******************************
            ARR2(I+1,NP)=ARR2(I,NP)
C******************************
11        CONTINUE
          I=0
2         ARR(I+1,NP)=A
C******************************
          ARR2(I+1,NP)=A2
C******************************
12      CONTINUE
        IF(JSTACK.EQ.0)RETURN
        IR=ISTACK(JSTACK)
        L=ISTACK(JSTACK-1)
        JSTACK=JSTACK-2
      ELSE
        K=(L+IR)/2
        TEMP=ARR(K,NP)
        ARR(K,NP)=ARR(L+1,NP)
        ARR(L+1,NP)=TEMP
C******************************
        TEMP2=ARR2(K,NP)
        ARR2(K,NP)=ARR2(L+1,NP)
        ARR2(L+1,NP)=TEMP2
C******************************
        IF(ARR(L+1,NP).GT.ARR(IR,NP))THEN
          TEMP=ARR(L+1,NP)
          ARR(L+1,NP)=ARR(IR,NP)
          ARR(IR,NP)=TEMP
C******************************
          TEMP2=ARR2(L+1,NP)
          ARR2(L+1,NP)=ARR2(IR,NP)
          ARR2(IR,NP)=TEMP2
C******************************
        ENDIF
        IF(ARR(L,NP).GT.ARR(IR,NP))THEN
          TEMP=ARR(L,NP)
          ARR(L,NP)=ARR(IR,NP)
          ARR(IR,NP)=TEMP
C******************************
          TEMP2=ARR2(L,NP)
          ARR2(L,NP)=ARR2(IR,NP)
          ARR2(IR,NP)=TEMP2
C******************************
        ENDIF
        IF(ARR(L+1,NP).GT.ARR(L,NP))THEN
          TEMP=ARR(L+1,NP)
          ARR(L+1,NP)=ARR(L,NP)
          ARR(L,NP)=TEMP
C******************************
          TEMP2=ARR2(L+1,NP)
          ARR2(L+1,NP)=ARR2(L,NP)
          ARR2(L,NP)=TEMP2
C******************************
        ENDIF
        I=L+1
        J=IR
        A=ARR(L,NP)
C******************************
        A2=ARR2(L,NP)
C******************************
3       CONTINUE
          I=I+1
        IF(ARR(I,NP).LT.A)GOTO 3
4       CONTINUE
          J=J-1
        IF(ARR(J,NP).GT.A)GOTO 4
        IF(J.LT.I)GOTO 5
        TEMP=ARR(I,NP)
        ARR(I,NP)=ARR(J,NP)
        ARR(J,NP)=TEMP
C******************************
        TEMP2=ARR2(I,NP)
        ARR2(I,NP)=ARR2(J,NP)
        ARR2(J,NP)=TEMP2
C******************************
        GOTO 3
5       ARR(L,NP)=ARR(J,NP)
        ARR(J,NP)=A
C******************************
        ARR2(L,NP)=ARR2(J,NP)
        ARR2(J,NP)=A2
C******************************
        JSTACK=JSTACK+2
        IF(JSTACK.GT.NSTACK) PRINT*,'NSTACK TOO SMALL IN SORT'
        IF(IR-I+1.GE.J-L)THEN
          ISTACK(JSTACK)=IR
          ISTACK(JSTACK-1)=I
          IR=J-1
        ELSE
          ISTACK(JSTACK)=J-1
          ISTACK(JSTACK-1)=L
          L=I
        ENDIF
      ENDIF
      GOTO 1
      END
