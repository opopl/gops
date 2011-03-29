!  GMIN: A PROGRAM FOR FINDING GLOBAL MINIMA
!  COPYRIGHT (C) 1999-2006 DAVID J. WALES
!  THIS FILE IS PART OF GMIN.
!
!  GMIN IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
!  IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
!  THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
!  (AT YOUR OPTION) ANY LATER VERSION.
!
!  GMIN IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
!  BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
!  GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
!
!  YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
!  ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
!  FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
!
!
!  THIS ROUTINE FINDS THE MATRIX R THAT MAPS A REFERENCE STRUCTURE
!  Q0 INTO A STRUCTURE Q, I.E.
!  Q = R Q0 FOR ALL N 3-DIMENSIONAL ATOMIC POSITION VECTORS OF Q AND Q0
!
SUBROUTINE REORIENT(NATOMS,Q,Q0,RMAT)
IMPLICIT NONE
INTEGER NATOMS, J1
DOUBLE PRECISION Q(3*NATOMS), Q0(3*NATOMS), RMAT(3,3), AMAT(3,3), CMAT(3,3), WORK(3), SERROR, RQ0(3)
DOUBLE PRECISION, PARAMETER :: EPS=1.0D-10
INTEGER IPIVOT(3), INFO

AMAT(1:3,1:3)=0.0D0
CMAT(1:3,1:3)=0.0D0
DO J1=1,NATOMS
   AMAT(1,1)=AMAT(1,1)+Q0(3*(J1-1)+1)*Q0(3*(J1-1)+1)+EPS ! EPS IS DESIGNED TO AVOID A SINGULAR MATRIX
   AMAT(1,2)=AMAT(1,2)+Q0(3*(J1-1)+1)*Q0(3*(J1-1)+2)+EPS ! FOR PLANAR SYSTEMS
   AMAT(1,3)=AMAT(1,3)+Q0(3*(J1-1)+1)*Q0(3*(J1-1)+3)+EPS
   AMAT(2,2)=AMAT(2,2)+Q0(3*(J1-1)+2)*Q0(3*(J1-1)+2)+EPS
   AMAT(2,3)=AMAT(2,3)+Q0(3*(J1-1)+2)*Q0(3*(J1-1)+3)+EPS
   AMAT(3,3)=AMAT(3,3)+Q0(3*(J1-1)+3)*Q0(3*(J1-1)+3)+EPS

   CMAT(1,1)=CMAT(1,1)+Q(3*(J1-1)+1)*Q0(3*(J1-1)+1)
   CMAT(1,2)=CMAT(1,2)+Q(3*(J1-1)+1)*Q0(3*(J1-1)+2)
   CMAT(1,3)=CMAT(1,3)+Q(3*(J1-1)+1)*Q0(3*(J1-1)+3)
   CMAT(2,1)=CMAT(2,1)+Q(3*(J1-1)+2)*Q0(3*(J1-1)+1)
   CMAT(2,2)=CMAT(2,2)+Q(3*(J1-1)+2)*Q0(3*(J1-1)+2)
   CMAT(2,3)=CMAT(2,3)+Q(3*(J1-1)+2)*Q0(3*(J1-1)+3)
   CMAT(3,1)=CMAT(3,1)+Q(3*(J1-1)+3)*Q0(3*(J1-1)+1)
   CMAT(3,2)=CMAT(3,2)+Q(3*(J1-1)+3)*Q0(3*(J1-1)+2)
   CMAT(3,3)=CMAT(3,3)+Q(3*(J1-1)+3)*Q0(3*(J1-1)+3)
ENDDO
! SYMMETRISE AMAT
AMAT(2,1)=AMAT(1,2); AMAT(3,1)=AMAT(1,3); AMAT(3,2)=AMAT(2,3)
! INVERT AMAT
CALL DGETRF(3,3,AMAT,3,IPIVOT,INFO) ! AMAT IS MODIFIED!
CALL DGETRI(3,AMAT,3,IPIVOT,WORK,3,INFO)
IF (INFO.NE.0) THEN
   PRINT '(A,I6)','ERROR - INFO AFTER DGETRI IN REORIENT=',INFO
   PRINT '(A)','Q:'
   PRINT '(3G20.10)',Q
   PRINT '(A)','Q0:'
   PRINT '(3G20.10)',Q0
   PRINT '(A)','AMAT:'
   PRINT '(3G20.10)',AMAT
   PRINT '(A)','CMAT:'
   PRINT '(3G20.10)',CMAT
   STOP
ENDIF

RMAT(1,1)=CMAT(1,1)*AMAT(1,1)+CMAT(1,2)*AMAT(2,1)+CMAT(1,3)*AMAT(3,1)
RMAT(1,2)=CMAT(1,1)*AMAT(1,2)+CMAT(1,2)*AMAT(2,2)+CMAT(1,3)*AMAT(3,2)
RMAT(1,3)=CMAT(1,1)*AMAT(1,3)+CMAT(1,2)*AMAT(2,3)+CMAT(1,3)*AMAT(3,3)
RMAT(2,1)=CMAT(2,1)*AMAT(1,1)+CMAT(2,2)*AMAT(2,1)+CMAT(2,3)*AMAT(3,1)
RMAT(2,2)=CMAT(2,1)*AMAT(1,2)+CMAT(2,2)*AMAT(2,2)+CMAT(2,3)*AMAT(3,2)
RMAT(2,3)=CMAT(2,1)*AMAT(1,3)+CMAT(2,2)*AMAT(2,3)+CMAT(2,3)*AMAT(3,3)
RMAT(3,1)=CMAT(3,1)*AMAT(1,1)+CMAT(3,2)*AMAT(2,1)+CMAT(3,3)*AMAT(3,1)
RMAT(3,2)=CMAT(3,1)*AMAT(1,2)+CMAT(3,2)*AMAT(2,2)+CMAT(3,3)*AMAT(3,2)
RMAT(3,3)=CMAT(3,1)*AMAT(1,3)+CMAT(3,2)*AMAT(2,3)+CMAT(3,3)*AMAT(3,3)

RETURN 
!
!  CHECK SQUARED ERROR - DEBUGGING ONLY.
!
SERROR=0.0D0
DO J1=1,NATOMS
   RQ0(1)=RMAT(1,1)*Q0(3*(J1-1)+1)+RMAT(1,2)*Q0(3*(J1-1)+2)+RMAT(1,3)*Q0(3*(J1-1)+3)
   RQ0(2)=RMAT(2,1)*Q0(3*(J1-1)+1)+RMAT(2,2)*Q0(3*(J1-1)+2)+RMAT(2,3)*Q0(3*(J1-1)+3)
   RQ0(3)=RMAT(3,1)*Q0(3*(J1-1)+1)+RMAT(3,2)*Q0(3*(J1-1)+2)+RMAT(3,3)*Q0(3*(J1-1)+3)
   SERROR=SERROR+(Q(3*(J1-1)+1)-RQ0(1))**2+(Q(3*(J1-1)+2)-RQ0(2))**2+(Q(3*(J1-1)+3)-RQ0(3))**2
ENDDO
PRINT '(A,G20.10)','ERROR IN REORIENT=',SERROR
IF (SERROR.GT.0.5D0) THEN
   RMAT(1:3,1:3)=-RMAT(1:3,1:3)
   SERROR=0.0D0
   DO J1=1,NATOMS
      RQ0(1)=RMAT(1,1)*Q0(3*(J1-1)+1)+RMAT(1,2)*Q0(3*(J1-1)+2)+RMAT(1,3)*Q0(3*(J1-1)+3)
      RQ0(2)=RMAT(2,1)*Q0(3*(J1-1)+1)+RMAT(2,2)*Q0(3*(J1-1)+2)+RMAT(2,3)*Q0(3*(J1-1)+3)
      RQ0(3)=RMAT(3,1)*Q0(3*(J1-1)+1)+RMAT(3,2)*Q0(3*(J1-1)+2)+RMAT(3,3)*Q0(3*(J1-1)+3)
      SERROR=SERROR+(Q(3*(J1-1)+1)-RQ0(1))**2+(Q(3*(J1-1)+2)-RQ0(2))**2+(Q(3*(J1-1)+3)-RQ0(3))**2
   ENDDO
   IF (SERROR.GT.0.5D0) THEN
      PRINT '(A,G20.10)','ERROR IN REORIENT FOR INVERSION=',SERROR
      STOP
   ENDIF
ENDIF

RETURN
END SUBROUTINE REORIENT
