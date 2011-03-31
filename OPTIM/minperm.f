C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
c     Interface to spjv.f for calculating minimum distance
c     of two atomic configurations with respect to
c     particle permutations.
c     The function permdist determines the distance or weight function,
c     minperm is the main routine.
c
c         Tomas Oppelstrup, Jul 10, 2003
c         tomaso@nada.kth.se
c

c     This is the main routine for minimum distance calculation.
c     Given two coordinate vectors p,q of particles each, return
c     the minimum distance in dist, and the permutation in perm.
c     perm is an integer vector such that
c       p(i) <--> q(perm(i))
c     i.e.
c       sum(i=1,n) permdist(p(i), q(perm(i))) == dist
c
      subroutine minperm(n, p, q, sx, sy, sz, pbc, perm, dist, worstdist, worstradius)
      implicit none

c     Input
c       n  : System size
c       p,q: Coordinate vectors (n particles)
c       s  : Box lengths (or dummy if open B.C.)
c       pbc: Periodic boundary conditions?
      integer n
      double precision p(3*n), q(3*n), s(3), sx, sy, sz, worstdist, worstradius
      logical pbc

c     Output
c       perm: Permutation so that p(i) <--> q(perm(i))
c       dist: Minimum attainable distance
c     We have
      integer perm(n)
      double precision dist, DUMMY
      
c     Parameters
c       scale : Precision
c       maxnei: Maximum number of closest neighbours
      double precision scale
      integer maxnei

      parameter (scale = 1.0d6   )
      parameter (maxnei = 60     )

c     Internal variables
c     cc, kk, first:
c       Sparse matrix of distances
c     first(i):
c       Beginning of row i in data,index vectors
c     kk(first(i)..first(i+1)-1):
c       Column indexes of existing elements in row i
c     cc(first(i)..first(i+1)-1):
c       Matrix elements of row i
      integer*8 kk(n*maxnei), first(n+1), x(n), y(n)
      integer*8 cc(n*maxnei), u(n), v(n), h
      integer   m, i, j, k, l, l2, t, a, i3, j3
      integer*8 n8, sz8, d
      INTEGER NDONE, J1, J2

c     Distance function
      double precision permdist

      s(1)=sx
      s(2)=sy
      s(3)=sz
      m = maxnei
      if(n .le. maxnei) m = n
      sz8 = m*n
      n8 = n

      do i=0,n
         first(i+1) = i*m + 1
      enddo
   
      if(m .eq. n) then
c     Compute the full matrix...

         do i=1,n
            k = first(i)-1
            do j=1,n
               cc(k+j) = permdist(p(3*i-2), q(3*j-2), s, pbc)*scale
               kk(k+j) = j
C              write(*,*) i, j, '-->', cc(k+j)
            enddo
         enddo
      else
c     We need to store the distances of the maxnei closeest neighbors
c     of each particle. The following builds a heap to keep track of
c     the maxnei closest neighbours seen so far. It might be more
c     efficient to use quick-select instead... (This is definately
c     true in the limit of infinite systems.)
        do i=1,n
           k = first(i)-1
           do j=1,m
              cc(k+j) = permdist(p(3*i-2), q(3*j-2), s, pbc)*scale
              kk(k+j) = j
              l = j
10            if(l .le. 1) goto 11
              l2 = l/2
              if(cc(k+l2) .lt. cc(k+l)) then
                 h = cc(k+l2)
                 cc(k+l2) = cc(k+l)
                 cc(k+l) = h
                 t = kk(k+l2)
                 kk(k+l2) = kk(k+l)
                 kk(k+l) = t
                 l = l2
                 goto 10
              endif
11         enddo
           
           do j=m+1,n
              d = permdist(p(3*i-2), q(3*j-2), s, pbc)*scale
              if(d .lt. cc(k+1)) then
                 cc(k+1) = d
                 kk(k+1) = j
                 l = 1
20               l2 = 2*l
                 if(l2+1 .gt. m) goto 21
                 if(cc(k+l2+1) .gt. cc(k+l2)) then
                    a = k+l2+1
                 else
                    a = k+l2
                 endif
                 if(cc(a) .gt. cc(k+l)) then
                    h = cc(a)
                    cc(a) = cc(k+l)
                    cc(k+l) = h
                    t = kk(a)
                    kk(a) = kk(k+l)
                    kk(k+l) = t
                    l = a-k
                    goto 20
                 endif
21               if (l2 .le. m) THEN ! split IF statements to avoid a segmentation fault
                    IF (cc(k+l2) .gt. cc(k+l)) then
                       h = cc(k+l2)
                       cc(k+l2) = cc(k+l)
                       cc(k+l) = h
                       t = kk(k+l2)
                       kk(k+l2) = kk(k+l)
                       kk(k+l) = t
                    ENDIF
                 endif
              endif
           enddo
!       PRINT '(A,I6,A)','atom ',i,' nearest neighbours and distances:'
!       PRINT '(20I6)',kk(m*(i-1)+1:m*i)
!       PRINT '(12I15)',cc(m*(i-1)+1:m*i)
           
        enddo    
!
! Create and maintain an ordered list, smallest to largest from kk(m*(i-1)+1:m*i) for atom i.
! NOTE that there is no symmetry with respect to exchange of I and J!
! This runs slower than the above heap algorithm.
!
!          CC(1:N*M)=HUGE(1)
!           DO I=1,N
!              K=FIRST(I)-1
!              DO J=1,N
!                 D=PERMDIST(P(3*I-2), Q(3*J-2), S, PBC)*SCALE
!                 IF (D.GT.CC(K+M)) CYCLE
!                 DO J1=M-1,1,-1
!                    IF (D.LT.CC(K+J1)) THEN
!                       CC(K+J1+1)=CC(K+J1)
!                       KK(K+J1+1)=KK(K+J1)
!                    ELSE
!                       CC(K+J1+1)=D
!                       KK(K+J1+1)=J
!                       GOTO 112
!                    ENDIF
!                 ENDDO
! !
! !  If we reach this point then we need to insert at the beginning.
! !
!                 CC(K+1)=D
!                 KK(K+1)=J
! 112             CONTINUE
!              ENDDO
!           ENDDO
      ENDIF

c     Call bipartite matching routine
      call jovosap(n8, sz8, cc, kk, first, x, y, u, v, h)

      if(h .lt. 0) then
c     If initial guess correct, deduce solution distance
c     which is not done in jovosap
         h = 0
         do i=1,n
            j = first(i)
 30         IF (J.GT.N*MAXNEI) THEN
!              PRINT '(A,I6,A)','minperm> WARNING A - matching failed'
               do J1=1,n
                  perm(J1)=J1
               enddo
               RETURN
            ENDIF
            if(kk(j) .ne. x(i)) then
               j = j + 1
               goto 30
            endif
            h = h + cc(j)
         enddo
      endif

      do i=1,n
         perm(i) = x(i)
         IF (PERM(I).GT.N) PERM(I)=N
         IF (PERM(I).LT.1) PERM(I)=1
      enddo

      dist = dble(h) / scale

      WORSTDIST=-1.0D0
      DO I=1,N
        DUMMY=(p(3*(i-1)+1)-q(3*(perm(i)-1)+1))**2+(p(3*(i-1)+2)-q(3*(perm(i)-1)+2))**2+(p(3*(i-1)+3)-q(3*(perm(i)-1)+3))**2
         IF (DUMMY.GT.WORSTDIST) THEN
            WORSTDIST=DUMMY 
            WORSTRADIUS=p(3*(i-1)+1)**2+p(3*(i-1)+2)**2+p(3*(i-1)+3)**2
         ENDIF
      ENDDO
      WORSTDIST=SQRT(WORSTDIST)
      WORSTRADIUS=MAX(SQRT(WORSTRADIUS),1.0D0)

      end
      
c     permdist is the distance or weight function. It is coded
c     separately for clarity. Just hope that the compiler
c     knows how to to do proper inlining!
c     Input
c       p,q: Coordinates
c       s  : Boxlengths (or dummy if open B.C.)
c       pbc: Periodic boundary conditions?

      double precision function permdist(p, q, s, pbc)
      implicit none
      double precision p(3), q(3), s(3)
      logical pbc

      double precision t, d
      integer i

      d = 0.0d0
      IF (PBC) THEN
         do i=1,3
            IF (S(I).NE.0.0D0) THEN
               t = q(i) - p(i)
               t = t - s(i)*anint(t/s(i))
               d = d + t*t
            ENDIF
         enddo
      ELSE
         d= (q(1) - p(1))**2+(q(2) - p(2))**2+(q(3) - p(3))**2
      ENDIF
      permdist = d
      end

c     The following routine performs weighted bipartite matching for
c     for a sparse non-negative integer weight matrix.
c     The original source is
c         http://www.magiclogic.com/assignment.html
c     A publication reference can be found on the above homepage and
c     in a comment below
c     

      SUBROUTINE JOVOSAP(N,SZ,CC,KK,FIRST,X,Y,U,V,H)
      IMPLICIT NONE
      INTEGER*8 N, SZ
      INTEGER*8 CC(SZ),KK(SZ),FIRST(N+1),X(N),Y(N),U(N),V(N)
      INTEGER*8 H,CNT,L0,T,T0,TD,V0,VJ,DJ
      INTEGER*8 LAB(N),D(N),FREE(N),TODO(N)
      LOGICAL OK(N)
      INTEGER*8 J, I, J0, L, J1, MIN, K, I0
      INTEGER*8 BIGINT

c     I don't know how to make g77 read integer*8 constants/parameters.
c       PARAMETER (BIGINT = 10**12) does not work(!)
c     nor does
c       PARAMETER (BIGINT = 1000000000000)
c     but this seems to be ok:
      BIGINT = 10**9
      BIGINT = BIGINT * 1000

C
C THIS SUBROUTINE SOLVES THE SPARSE LINEAR ASSIGNMENT PROBLEM
C ACCORDING 
C
C   "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear   
C    Assignment Problems," Computing 38, 325-340, 1987
C   
C   by
C   
C   R. Jonker and A. Volgenant, University of Amsterdam.
C
C
C INPUT PARAMETERS :
C N = NUMBER OF ROWS AND COLUMNS
C C = WEIGHT MATRIX
C
C OUTPUT PARAMETERS
C X = COL ASSIGNED TO ROW
C Y = ROW ASSIGNED TO COL
C U = DUAL ROW VARIABLE
C V = DUAL COLUMN VARIABLE
C H = VALUE OF OPTIMAL SOLUTION
C
C INITIALIZATION

c     Next line added by tomaso@nada.kth.se, to enable detection
c     of solutions being equivalent to the initial guess

C
C  If Y(:) is initialised to zero then we see segmentation faults if 
C  a Y element is unset, etc.
C

      Y(1:N) = 0
      X(1:N) = 0
      TODO(1:N)=0
      h = -1
      DO 10 J=1,N
         V(J)=BIGINT
   10 CONTINUE
      DO 20 I=1,N
         X(I)=0
         DO 15 T=FIRST(I),FIRST(I+1)-1
            J=KK(T)
            IF (CC(T).LT.V(J)) THEN
              V(J)=CC(T)
              Y(J)=I
            END IF
   15    CONTINUE
   20 CONTINUE
      DO 30 J=1,N
         J0=N-J+1
         I=Y(J0)
         IF (I.EQ.0) THEN
!           PRINT '(A,I6,A)','minperm> WARNING B - matching failed'
            RETURN
         ENDIF
         IF (X(I).NE.0) THEN
           X(I)=-ABS(X(I))
           Y(J0)=0
         ELSE
           X(I)=J0
         END IF
   30 CONTINUE
      L=0
      DO 40 I=1,N
         IF (X(I).EQ.0) THEN
           L=L+1
           FREE(L)=I
           GOTO 40
         END IF
         IF (X(I).LT.0) THEN
           X(I)=-X(I)
         ELSE
           J1=X(I)
           MIN=BIGINT
           DO 31 T=FIRST(I),FIRST(I+1)-1
              J=KK(T)
              IF (J.EQ.J1) GOTO 31
              IF (CC(T)-V(J).LT.MIN) MIN=CC(T)-V(J)
   31      CONTINUE
           V(J1)=V(J1)-MIN
         END IF
   40 CONTINUE
C IMPROVE THE INITIAL SOLUTION
      CNT=0
      IF (L.EQ.0) GOTO 1000
   41 L0=L
      K=1
      L=0
   50 I=FREE(K)
      K=K+1
      V0=BIGINT
      VJ=BIGINT
      DO 42 T=FIRST(I),FIRST(I+1)-1
         J=KK(T)
         H=CC(T)-V(J)
         IF (H.LT.VJ) THEN
           IF (H.GE.V0) THEN
             VJ=H
             J1=J
           ELSE
             VJ=V0
             V0=H
             J1=J0
             J0=J
           END IF
         END IF
   42 CONTINUE
      I0=Y(J0)
      IF (V0.LT.VJ) THEN
        V(J0)=V(J0)-VJ+V0
      ELSE
        IF (I0.EQ.0) GOTO 43
        J0=J1
        I0=Y(J1)
      END IF
      IF (I0.EQ.0) GOTO 43
      IF (V0.LT.VJ) THEN
        K=K-1
        FREE(K)=I0
      ELSE
        L=L+1
        FREE(L)=I0
      END IF
   43 X(I)=J0
      Y(J0)=I
      IF (K.LE.L0) GOTO 50
      CNT=CNT+1
      IF ((L.GT.0).AND.(CNT.LT.2)) GOTO 41
C AUGMENTATION PART
      L0=L
      DO 90 L=1,L0
         DO 51 J=1,N
            OK(J)=.FALSE.
            D(J)=BIGINT
   51    CONTINUE
         MIN=BIGINT
         I0=FREE(L)
         TD=N
         DO 52 T=FIRST(I0),FIRST(I0+1)-1
            J=KK(T)
            DJ=CC(T)-V(J)
            D(J)=DJ
            LAB(J)=I0
            IF (DJ.LE.MIN) THEN
              IF (DJ.LT.MIN) THEN
                MIN=DJ
                K=1
                TODO(1)=J
              ELSE
                K=K+1
                TODO(K)=J
              END IF
            END IF
   52    CONTINUE
         DO 53 H=1,K
            J=TODO(H)
            IF (J.EQ.0) THEN
!              PRINT '(A,I6,A)','minperm> WARNING C - matching failed'
               RETURN
            ENDIF
            IF (Y(J).EQ.0) GOTO 80
            OK(J)=.TRUE.
   53    CONTINUE
C REPEAT UNTIL A FREE ROW HAS BEEN FOUND
   60    IF (K.EQ.0) THEN
!           PRINT '(A,I6,A)','minperm> WARNING D - matching failed'
            RETURN
         ENDIF
         J0=TODO(K)
         K=K-1
         I=Y(J0)
         TODO(TD)=J0
         TD=TD-1
         T0=FIRST(I)
         T=T0-1
   61    T=T+1
         IF (KK(T).NE.J0) GOTO 61
         H=CC(T)-V(J0)-MIN
         DO 62 T=T0,FIRST(I+1)-1
            J=KK(T)
            IF (.NOT. OK(J)) THEN
              VJ=CC(T)-H-V(J)
              IF (VJ.LT.D(J)) THEN
                D(J)=VJ
                LAB(J)=I
                IF (VJ.EQ.MIN) THEN
                  IF (Y(J).EQ.0) GOTO 70
                  K=K+1
                  TODO(K)=J
                  OK(J)=.TRUE.
                END IF
              END IF
            END IF
   62    CONTINUE
         IF (K.NE.0) GOTO 60
         MIN=BIGINT-1
         DO 63 J=1,N
            IF (D(J).LE.MIN) THEN
              IF (.NOT. OK(J)) THEN
                IF (D(J).LT.MIN) THEN
                  MIN=D(J)
                  K=1
                  TODO(1)=J
                ELSE
                  K=K+1
                  TODO(K)=J
                END IF
              END IF
            END IF
   63    CONTINUE
         DO 64 J0=1,K
            J=TODO(J0)
            IF (Y(J).EQ.0) GOTO 70
            OK(J)=.TRUE.
   64    CONTINUE
         GOTO 60
   70    IF (MIN.EQ.0) GOTO 80
         DO 71 K=TD+1,N
            J0=TODO(K)
            V(J0)=V(J0)+D(J0)-MIN
   71    CONTINUE
   80    I=LAB(J)
         Y(J)=I
         K=J
         J=X(I)
         X(I)=K
         IF (I0.NE.I) GOTO 80
   90 CONTINUE
      H=0
      DO 100 I=1,N
         J=X(I)
         T=FIRST(I)-1
  101    T=T+1
         IF (T.GT.SZ) THEN
            PRINT '(A,I6,A)','minperm> WARNING D - atom ',I,' not matched - maximum number of neighbours too small?'
            RETURN
         ENDIF
         IF (KK(T).NE.J) GOTO 101
         DJ=CC(T)
         U(I)=DJ-V(J)
         H=H+DJ
  100 CONTINUE

 1000 END
