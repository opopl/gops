!op226> GPL License Info {{{
!  GMIN: A program for finding global minima
!  Copyright (C) 1999-2006 David J. Wales
!  This file is part of GMIN.
!
!  GMIN is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  GMIN is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!op226>}}}

MODULE NOA
      IMPLICIT NONE
      SAVE
      INTEGER :: NUMBER_OF_ATOMS

      CONTAINS

      SUBROUTINE COUNTATOMS(MYUNIT)
!op226> Declarations {{{ 
      IMPLICIT NONE
      INTEGER :: EOF,NRES,SEQ(500),I_RES,NOGLY,GLY, MYUNIT
      LOGICAL :: YESNO
      CHARACTER(LEN=5) TARFL
      CHARACTER(LEN=10)  CHECK
      CHARACTER(LEN=80) MYLINE,INPCRD1
!op226>}}} 

      YESNO=.FALSE.

      INQUIRE(FILE='coords',EXIST=YESNO)

      NUMBER_OF_ATOMS=0
      ! read in coordinate files:
      ! coords  
      ! pro.list        amh
      IF (YESNO) THEN
         OPEN(UNIT=7,FILE='coords',STATUS='OLD')
         DO
            READ(7,*,IOSTAT=EOF)
            IF (EOF==0) THEN
               NUMBER_OF_ATOMS = NUMBER_OF_ATOMS + 1
            ELSE
               EXIT
            ENDIF
         ENDDO
            ELSE
         PRINT '(A)','ERROR - no coords file'
         STOP
      ENDIF
      CLOSE(7)
      
      END SUBROUTINE COUNTATOMS

END MODULE NOA
