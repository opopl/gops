	SUBROUTINE  GENTAB_RAMA

        USE GLOBALS,ONLY:AMINOA,MAXSIZ,TGSEQUENCES_AMC,RAMA_PROB, &
                 RAMASCL,RAMA_FORCE

	IMPLICIT NONE

!     INTERNAL VARIABLES:

!   RAMA_PROB(PHI_INDEX,PSI_INDEX,AMINO,AMINO-1,AMINO,AMINO+1

        INTEGER ISIT1,ISIT2,MVM_SCR
        INTEGER IAA,IRES,IPRE,IPOST,OPEN_STATUS,I1,NMRES

!     REQUIRED SUBROUTINES

!       EXTERNAL 

!  DATA IN 10 DEGREE WIDE BINS FROM 180 TO -180
 
        RAMA_PROB(:,:,:,:,:)=0.0

!!$        DO IAA = 1,20
!!$         DO IPRE = 1,20
!!$          DO IPOST = 1,20
!!$
!!$        OPEN(MVM_SCR, &
!!$          FILE='~/MARCIO/BDTRIMERS/'//AMINOA(IAA)//'/'//AMINOA(IPRE)//'_'//AMINOA(IAA)//'_'//AMINOA(IPOST)//'.SCR',STATUS='OLD',IOSTAT=OPEN_STATUS)
!!$               IF (OPEN_STATUS.NE.0) THEN
!!$                 WRITE(6,*) 'FAILURE TO OPEN FILE IN GENTAB_RAMA'
!!$                 STOP
!!$               ENDIF
!!$        DO ISIT1 = 1,36
!!$        READ(MVM_SCR,*)(RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST),ISIT2=1,36)
!!$        ENDDO
!!$        CLOSE(MVM_SCR)
!!$ 
!!$          ENDDO ! DO IPOST = 1,20
!!$         ENDDO ! DO IPRE = 1,20
!!$        ENDDO ! DO IAA = 1,20
!!$
!!$!   LOCAL FILTERING  
!!$
!!$        DO IAA = 1,20
!!$         DO IPRE = 1,20
!!$          DO IPOST = 1,20
!!$             DO ISIT1 = 1,36
!!$              DO ISIT2 = 1,36
!!$
!!$               IF(RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST).GT.0.0)THEN
!!$                RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = & 
!!$                -RAMASCL* LOG(RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST))
!!$               ENDIF 
!!$
!!$                IF((ISIT1.EQ.1).AND.(ISIT2.EQ.1))THEN
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2,IPRE,IAA,IPOST) 
!!$
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) / 4
!!$               ELSE IF((ISIT1.EQ.1).AND.(ISIT2.EQ.36))THEN
!!$                 RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1,ISIT2-1,IPRE,IAA,IPOST)
!!$
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) / 4
!!$               ELSE IF((ISIT1.EQ.36).AND.(ISIT2.EQ.36))THEN
!!$
!!$                 RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2,IPRE,IAA,IPOST) 
!!$
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) / 4
!!$
!!$               ELSE IF((ISIT1.EQ.36).AND.(ISIT2.EQ.1))THEN
!!$                 RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2,IPRE,IAA,IPOST)  
!!$
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST)/ 4
!!$               ELSE IF((ISIT1.EQ.1).AND.((ISIT2.NE.1).OR.(ISIT2.NE.36)))THEN
!!$                 RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2,IPRE,IAA,IPOST)
!!$
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) / 6
!!$          ELSE IF((ISIT1.EQ.36).AND.((ISIT2.NE.1).OR.(ISIT2.NE.36)))THEN
!!$                 RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2,IPRE,IAA,IPOST) 
!!$
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) / 6
!!$          ELSE IF((ISIT2.EQ.1) .AND. ((ISIT1.NE.1).OR.(ISIT1.NE.36)))THEN
!!$                 RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2,IPRE,IAA,IPOST)
!!$
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) / 6
!!$          ELSE IF((ISIT2.EQ.36).AND.((ISIT1.NE.1).OR.(ISIT1.NE.36)))THEN
!!$                 RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2,IPRE,IAA,IPOST)
!!$
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) / 6
!!$               ELSE 
!!$                 RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1,ISIT2+1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1,ISIT2-1,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1-1,ISIT2,IPRE,IAA,IPOST) +  &
!!$                  RAMA_PROB(ISIT1+1,ISIT2,IPRE,IAA,IPOST)
!!$
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) / 9
!!$               ENDIF
!!$
!!$              ENDDO  ! DO ISIT1 = 1,36
!!$            ENDDO ! DO ISIT2 = 1,36
!!$          ENDDO ! DO IPOST = 1,20
!!$         ENDDO ! DO IPRE = 1,20
!!$        ENDDO ! DO IAA = 1,20
!!$
!!$! CALC RAMA_FORCE
!!$
!!$          DO IAA = 1,20
!!$           DO IPRE = 1,20
!!$            DO IPOST = 1,20
!!$             DO ISIT1 = 1,36
!!$              DO ISIT2 = 1,36
!!$
!!$             IF ( ISIT2 .EQ. 1) THEN
!!$               RAMA_FORCE(ISIT1,ISIT2,IPRE,IAA,IPOST) =  & 
!!$                  (RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) - & 
!!$                   RAMA_PROB(ISIT1,36,IPRE,IAA,IPOST) )/2
!!$             ELSE IF ( ISIT2 .EQ. 1) THEN 
!!$               RAMA_FORCE(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                  (RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) - &
!!$                   RAMA_PROB(36,ISIT2,IPRE,IAA,IPOST) )/2
!!$             ELSE 
!!$               RAMA_FORCE(ISIT1,ISIT2,IPRE,IAA,IPOST) = &
!!$                 (RAMA_PROB(ISIT1,ISIT2,IPRE,IAA,IPOST) - & 
!!$                  RAMA_PROB(ISIT1,ISIT2-1,IPRE,IAA,IPOST) )/2
!!$             ENDIF
!!$
!!$              ENDDO ! DO ISIT2 = 1,36
!!$             ENDDO ! DO ISIT1 = 1,36
!!$            ENDDO ! DO IPOST = 1,20
!!$           ENDDO ! DO IPRE = 1,20
!!$          ENDDO ! DO IAA = 1,20

        RETURN
	END
