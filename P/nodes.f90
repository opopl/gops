!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

MODULE NODES
     USE COMMONS,ONLY : SLURMT
     IMPLICIT NONE
     SAVE

     INTEGER :: JPN,NNODES,NNLEN
     CHARACTER(LEN=80),ALLOCATABLE,DIMENSION(:) :: NODENAME
     CHARACTER(LEN=80) :: USERNAME
     CHARACTER(LEN=100) :: WORKINGDIRECTORY
     CHARACTER(LEN=1000) :: TOTALJOBSTRING
     CHARACTER(LEN=200) :: PATHSTRING
     LOGICAL :: SSHPARALLEL=.FALSE., YESNO

     CHARACTER(LEN=120) :: NODESTRING
     CHARACTER(LEN=120) :: NODEN
     CHARACTER(LEN=120) :: TEMPSTRING
     INTEGER NSTART, NFINISH, NSTART2, NFINISH2, J1, N1, N2, LCOUNT

     CONTAINS

     SUBROUTINE GETNODES(NCPU)
          USE PORFUNCS
          IMPLICIT NONE
     
          INTEGER,INTENT(OUT) :: NCPU

          INTEGER :: I
          INTEGER,PARAMETER :: UNIT1=91
          CHARACTER(LEN=80) :: ARG

          WRITE(*,'(a,i1,a)') 'GetNodes> Will run ',JPN,' jobs per node'
          INQUIRE(FILE='nodes.info',EXIST=YESNO)
          IF (YESNO) THEN
             OPEN(UNIT=1,FILE='nodes.info',STATUS='OLD')
             IF (SLURMT) THEN
                READ(1,'(A)') NODESTRING
                READ(1,*) NNODES
                IF (NNODES.EQ.1) THEN ! no list to pick apart!
                   NNLEN=JPN*NNODES
                   IF (ALLOCATED(NODENAME)) DEALLOCATE(NODENAME) ! to allow for calling keyword more than once
                   ALLOCATE(NODENAME(NNLEN))
                   NCPU=JPN*NNODES
                   WRITE(*,'(A,I2)') 'GetNodes> NCPU = ',NCPU
                   SSHPARALLEL=.TRUE.   
                   WRITE(*,'(A,I2,A)') 'GetNodes> Following ',NNODES,' nodes are available:'
                   PRINT '(A)', TRIM(ADJUSTL(NODESTRING))
                   NODENAME(1:NCPU)= TRIM(ADJUSTL(NODESTRING))
                ELSE
                   LCOUNT=0
                   DO 
                      LCOUNT=LCOUNT+1
                      IF (NODESTRING(LCOUNT:LCOUNT).EQ.'[') THEN
                         NODEN=NODESTRING(1:LCOUNT-1)
                         EXIT
                      ENDIF
                   ENDDO 
                   NNLEN=JPN*NNODES
                   IF (ALLOCATED(NODENAME)) DEALLOCATE(NODENAME) ! to allow for calling keyword more than once
                   ALLOCATE(NODENAME(NNLEN))
                   NCPU=JPN*NNODES
                   WRITE(*,'(A,I2)') 'GetNodes> NCPU = ',NCPU
                   SSHPARALLEL=.TRUE.   
                   WRITE(*,'(A,I2,A)') 'GetNodes> Following ',NNODES,' nodes are available:'
                   NSTART=LCOUNT+1
                   I=0
                   outer: DO 
                      LCOUNT=LCOUNT+1
                      IF (NODESTRING(LCOUNT:LCOUNT).EQ.'-') THEN
                         NFINISH=LCOUNT-1
                         NSTART2=LCOUNT+1
                         DO 
                            LCOUNT=LCOUNT+1
                            IF ((NODESTRING(LCOUNT:LCOUNT).EQ.',').OR.(NODESTRING(LCOUNT:LCOUNT).EQ.']')) THEN
                               NFINISH2=LCOUNT-1
   
                               READ(NODESTRING(NSTART:NFINISH),*) N1
                               READ(NODESTRING(NSTART2:NFINISH2),*) N2
   !                           PRINT '(A,6I8)','NSTART,NFINISH,NSTART2,NFINISH2,N1,N2=',NSTART,NFINISH,NSTART2,NFINISH2,N1,N2
                               
                               DO J1=N1,N2
                                  I=I+1
                                  WRITE(TEMPSTRING,'(I10)') J1
                                  NODENAME(JPN*(I-1)+1:JPN*I)= TRIM(ADJUSTL(NODEN)) // TRIM(ADJUSTL(TEMPSTRING))
                                  PRINT '(A)', TRIM(ADJUSTL(NODEN))  // TRIM(ADJUSTL(TEMPSTRING))
                               ENDDO
   
                               NSTART=LCOUNT+1
                               IF (NODESTRING(LCOUNT:LCOUNT).EQ.']') EXIT outer
                               EXIT 
                            ENDIF
                         ENDDO
                      ELSEIF (NODESTRING(LCOUNT:LCOUNT).EQ.',') THEN
                          NFINISH=LCOUNT-1
                          I=I+1
                          NODENAME(JPN*(I-1)+1:JPN*I)= TRIM(ADJUSTL(NODEN)) // NODESTRING(NSTART:NFINISH)
                          PRINT '(A)', TRIM(ADJUSTL(NODEN)) // NODESTRING(NSTART:NFINISH)
                          NSTART=LCOUNT+1
                      ELSEIF (NODESTRING(LCOUNT:LCOUNT).EQ.']') THEN
                          NFINISH=LCOUNT-1
                          I=I+1
                          NODENAME(JPN*(I-1)+1:JPN*I)= TRIM(ADJUSTL(NODEN)) // NODESTRING(NSTART:NFINISH)
                          PRINT '(A)', TRIM(ADJUSTL(NODEN)) // NODESTRING(NSTART:NFINISH)
                         EXIT outer
                      ENDIF
                   ENDDO outer
                ENDIF
                READ(1,'(A)') USERNAME
                READ(1,'(A)') WORKINGDIRECTORY
                WRITE(*,'(2A)') 'GetNodes> Working in directory ',TRIM(ADJUSTL(WORKINGDIRECTORY))
                CLOSE(1)
                PRINT '(A,I8,A)','GetNodes> Complete list of ',NCPU,' cores:'
                PRINT '(A)',NODENAME(1:NCPU)
             ELSE
                READ(1,*) NNODES
                NNLEN=JPN*NNODES
                IF (ALLOCATED(NODENAME)) DEALLOCATE(NODENAME) ! to allow for calling keyword more than once
                ALLOCATE(NODENAME(NNLEN))
                NCPU=JPN*NNODES
                WRITE(*,'(a,i2)') 'GetNodes> NCPU = ',ncpu
             
                SSHPARALLEL=.TRUE.   
                WRITE(*,'(a,i2,a)') 'GetNodes> Following ',NNODES,' nodes are available:'
                DO I=1,NNODES
                     READ(1,'(a)') arg
                     PRINT '(A)', TRIM(ADJUSTL(arg))
                     NODENAME(JPN*(I-1)+1:JPN*I)=ADJUSTL(ARG)
                ENDDO
                READ(1,'(A)') USERNAME
                READ(1,'(A)') WORKINGDIRECTORY
                WRITE(*,'(2A)') 'GetNodes> Working in directory ',TRIM(ADJUSTL(WORKINGDIRECTORY))
                CLOSE(1)
             ENDIF
          ELSE
             PRINT '(A)','getnodes> No nodes.info file - assuming no OPTIM jobs required or interactive run'
          ENDIF
     END SUBROUTINE GETNODES

     SUBROUTINE SSHSUBMIT(ICPU,STAT,JOBSTRING,CONNSTR1,LDEBUG)
          USE PORFUNCS, ONLY: SYSTEM_SUBR
          USE COMMONS, ONLY: CHARMMT, ZSYM, COPYFILES, COPYOPTIMT, BHINTERPT, BISECTT
          IMPLICIT NONE

          INTEGER,INTENT(IN) :: ICPU
          LOGICAL,INTENT(IN) :: LDEBUG
          INTEGER,INTENT(OUT) :: STAT
          CHARACTER(LEN=*) :: JOBSTRING,CONNSTR1

          INTEGER :: MYSTATUS
          CHARACTER(LEN=80) :: NODE

          NODE=TRIM(ADJUSTL(NODENAME(ICPU)))
!
!  There is a problem here with the return status. We cannot get the return status
!  from the OPTIM job - instead we get the return status from rsh, which will be 0 unless
!  something in the job submission actually fails. If we have a path.info.connstr1 file
!  then we can analyse it, so this rsh status is actually the important one!
!
!  Job submission now changed to use a single system call and one large job string.
!  Putting all the data transfer etc. in the original rsh, so that it runs on
!  the compute node, should avoid any other rcp, rsync or rsh !
!
          IF (SLURMT) THEN

             PATHSTRING=TRIM(ADJUSTL(WORKINGDIRECTORY)) // '/' // connstr1
!
! Build up the complete rsh command step by step:
! (1) make the scratch directory on the node. -p flag means no error is generated if the directory already exists.
             TOTALJOBSTRING= 'ssh ' // TRIM(node) // ' " mkdir -p ' // TRIM(ADJUSTL(PATHSTRING)) 
! (2) move to the WORKINGDIRECTORY (saves unpicking the COPYFILES list!)
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' && cd ' // TRIM(ADJUSTL(WORKINGDIRECTORY)) 
! (3) copy data from WORKINGDIRECTORY to the scratch directory on the node. 
!     Note that if any file is missing an error condition will result, and subsequent commands will fail.
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
  &             // ' && cp -r  *.' // connstr1 // ' ' // TRIM(ADJUSTL(COPYFILES)) // ' ' // TRIM(ADJUSTL(PATHSTRING))
! (4) move to the scratch directory on the node
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' && cd ' // TRIM(ADJUSTL(PATHSTRING))
! (4b) delete any existing path.info.* file (a very rare but not impossible condition!)
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' && rm -f path.info.* '
! (5) run the OPTIM job
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' && ' // JOBSTRING
! (6) copy results back
             IF (LDEBUG) THEN ! copy everything back
                TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' && cp  *.' // connstr1 &
   &                       // ' ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
             ELSEIF (COPYOPTIMT.AND.(BHINTERPT.OR.BISECTT)) THEN ! copy path.info, OPTIM
                TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &               // ' && cp OPTIM* min.data.info* ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
             ELSEIF (COPYOPTIMT) THEN ! copy path.info, OPTIM, odata and finish
                TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &               // ' && cp OPTIM* path.info* ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
             ELSE ! we only really need path.info or min.data.info
                IF (BHINTERPT.OR.BISECTT) THEN
                   TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &               // ' && cp min.data.info* ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
                ELSE
                   TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &               // ' && cp path.info* ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
                ENDIF
             ENDIF
! (7) remove the scratch directory
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' && rm -r ' // TRIM(ADJUSTL(PATHSTRING)) // ' " '
             IF (LDEBUG) PRINT '(2A)', 'nodes> complete job string: ',TRIM(ADJUSTL(TOTALJOBSTRING))
! (8) submit the job for real
             CALL SYSTEM_SUBR(TRIM(ADJUSTL(TOTALJOBSTRING)),MYSTATUS)  
             STAT=MYSTATUS

          ELSE ! rsh works on the other clusters

             PATHSTRING='/scratch/' // TRIM(ADJUSTL(USERNAME)) // '/' // connstr1
!
! Build up the complete rsh command step by step:
! (1) make the scratch directory on the node. -p flag means no error is generated if the directory already exists.
             TOTALJOBSTRING= 'rsh ' // TRIM(node) // ' " mkdir -p ' // TRIM(ADJUSTL(PATHSTRING)) 
! (2) move to the WORKINGDIRECTORY (saves unpicking the COPYFILES list!)
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; cd ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
! (3) copy data from WORKINGDIRECTORY to the scratch directory on the node. 
!     Note that if any file is missing an error condition will result, and subsequent commands will fail.
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
  &             // ' ; cp -r  *.' // connstr1 // ' ' // TRIM(ADJUSTL(COPYFILES)) // ' ' // TRIM(ADJUSTL(PATHSTRING))
! (4) move to the scratch directory on the node
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; cd ' // TRIM(ADJUSTL(PATHSTRING))
! (4b) delete any existing path.info.* file (a very rare but not impossible condition!)
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' && rm -f path.info.* '
! (5) run the OPTIM job
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; ' // JOBSTRING
! (6) copy results back
             IF (LDEBUG) THEN ! copy everything back 
                TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; cp  *.' // connstr1 &
   &                         // ' ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
             ELSEIF (COPYOPTIMT.AND.(BHINTERPT.OR.BISECTT)) THEN ! copy path.info, OPTIM
                TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &               // ' ; cp OPTIM* min.data.info* ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
             ELSEIF (COPYOPTIMT) THEN ! copy path.info, OPTIM, odata and finish
                TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &               // ' ; cp OPTIM* path.info* ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
             ELSE ! we only really need path.info or min.data.info
                IF (BHINTERPT.OR.BISECTT) THEN
                   TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &               // ' ; cp min.data.info* ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
                ELSE
                   TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
   &               // ' ; cp path.info* ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
                ENDIF
             ENDIF
! (7) remove the scratch directory
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; rm -r ' // TRIM(ADJUSTL(PATHSTRING)) // ' " '
             IF (LDEBUG) PRINT '(2A)', 'nodes> complete job string: ',TRIM(ADJUSTL(TOTALJOBSTRING)) 
! (8) submit the job for real
             CALL SYSTEM_SUBR(TRIM(ADJUSTL(TOTALJOBSTRING)),MYSTATUS)  
             STAT=MYSTATUS
          ENDIF
     END SUBROUTINE SSHSUBMIT
END MODULE NODES
