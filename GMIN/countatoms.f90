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
      USE MODAMBER9 , ONLY : INPCRD
      IMPLICIT NONE
      SAVE
      INTEGER :: NUMBER_OF_ATOMS

      CONTAINS

      SUBROUTINE COUNTATOMS(MYUNIT)
!op226> Declarations {{{ 
      IMPLICIT NONE
      INTEGER :: EOF,NRES,SEQ(500),I_RES,NOGLY,GLY, MYUNIT
      LOGICAL :: YESNO, YESNOA, YESNOC, YESNOAMH, YESNOA9
      CHARACTER(LEN=5) TARFL
      CHARACTER(LEN=10)  CHECK
      CHARACTER(LEN=80) MYLINE,INPCRD1
!op226>}}} 

      YESNO=.FALSE.
      YESNOA=.FALSE.
      YESNOAMH=.FALSE.
      YESNOA9=.FALSE.
!
!  If the current working directory contains more than one of these files
!  then the precedence is coords, then input.crd, then coords.amber
!  OPTIM does this a bit better by calling main.f first to see if
!  we are actually doing AMBER or CHARMM. 
!
      INQUIRE(FILE='pro.list',EXIST=YESNOAMH)
      INQUIRE(FILE='coords',EXIST=YESNO)
      INQUIRE(FILE='coords.amber',EXIST=YESNOA)
      INQUIRE(FILE='input.crd',EXIST=YESNOC)
      INQUIRE(FILE='coords.inpcrd',EXIST=YESNOA9)

      NUMBER_OF_ATOMS=0
      ! read in coordinate files:
      ! coords  
      ! pro.list        amh
      ! input.crd       charmm
      ! coords.inpcrd   amber
      ! coords.amber    
      !     {{{
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
        ELSEIF (YESNOAMH) THEN
         open(unit=30,file='pro.list',status='old',form='formatted')
         read (30,1000)tarfl
1000     format(a5)
         close(30)

          open(30,file='proteins/'//tarfl,status='old')
            read(30,*)
            read(30,*)nres
            if (nres.gt.500) then
                write(6,*) 'failure nres gr than 500 countatoms'
                stop
            endif
            read (30,25)(seq(i_res),i_res=1,nres)
!            write(6,25)(seq(i_res),i_res=1,nres)
25         format(25(i2,1x))
          close(30)

          NOGLY = 0
          GLY = 0

           do i_res=1,nres
             if (seq(i_res).ne.8) NOGLY = NOGLY +1
             if (seq(i_res).eq.8) GLY = GLY +1
           enddo

            Number_of_Atoms = NOGLY*3 + GLY*2
      ELSE IF (YESNOA9) THEN
!         OPEN(UNIT=7,FILE='coords.gayberne',STATUS='OLD')
!         PRINT '(A)','reading coordinates from file coords.gayberne'

         inpcrd1='coords.inpcrd'
!         inpcrd1=trim(adjustl(inpcrd1))
         call amberinterface(Number_of_Atoms,1,inpcrd1,MYUNIT)

      ELSEIF (YESNOC) THEN
         OPEN(UNIT=7,FILE='input.crd',STATUS='OLD')
         do
           read(7,*) myline
           if (myline(1:1)=='*') then ! SAT This is the goddamn CHARMM comment line
              cycle
           else
              read(myline,*) Number_of_Atoms
              exit
           endif
         enddo

! DAE We also need to find out what MAXAIM is in CHARMM, and set MXATMS in OPTIM to be the same, so that those arrays which
! are passed between the two can be declared correctly. MXATMS is now stored in modmxatms.

         CALL GETMAXAIM
         WRITE(MYUNIT,'(A,I8)') 'countatoms> Number_of_Atoms=',Number_of_Atoms
      ELSEIF (YESNOA) THEN
         OPEN(UNIT=7,FILE='coords.amber',STATUS='OLD')
         do
            read(7,'(A3)',iostat=eof) check
            if (eof.LT.0) then
               PRINT *,'End of file before all information specified'
               STOP
            ENDIF
            IF (check.EQ.'end' .OR. check.EQ.'END' .OR. check.EQ.'End') THEN
               CLOSE(7)
               EXIT
            ENDIF
            Number_of_Atoms = Number_of_Atoms + 1
         enddo
      ELSE
         PRINT '(A)','ERROR - no coords, input.crd, coords.inpcrd or coords.amber file'
         STOP
      ENDIF
      ! }}}
      CLOSE(7)

!     print *, Number_of_Atoms, ' atoms in the system'
      
      END SUBROUTINE COUNTATOMS

END MODULE NOA
