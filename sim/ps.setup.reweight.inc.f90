!
!     Optional change of reactant minima set via reweighting.
!
      IF (REWEIGHTT) THEN
        ! {{{
         ALLOCATE(CANDIDATES(NMIN))
         IF (DIRECTION.EQ.'AB') THEN
            ALLOCATE(NEWPFMIN(NMINB))
            NEWPFMIN(1:NMINB)=0.0D0
!
!  Select NRWREACTANT minima from the B set according to the required weights in RWPROB
!
            PFTOTALB=0.0D0
            DO J1=1,NRWBINS
               IF (NINT(NRWREACTANT*RWPROB(J1)).EQ.0) CYCLE
               NCOUNT=0
               DO J2=1,NMINB ! identify minima in the required energy range
                  IF ((EMIN(LOCATIONB(J2)).GE.RWEMIN+(J1-1)*RWBINWIDTH).AND.
     &                (EMIN(LOCATIONB(J2)).LE.RWEMIN+J1*RWBINWIDTH)) THEN
                      NCOUNT=NCOUNT+1
                      CANDIDATES(NCOUNT)=J2
                  ENDIF
               ENDDO
               PRINT '(3(A,I8),A,G20.10)','setup> number of B minima in energy bin ',J1,' is ',NCOUNT,' number needed=',
     &                           NINT(NRWREACTANT*RWPROB(J1)),' probability=',RWPROB(J1)
               IF (NCOUNT.EQ.0) STOP
               DO J2=1,NINT(NRWREACTANT*RWPROB(J1))
                  NRANDOM=NCOUNT*DPRAND()+1
                  PRINT '(3(A,I8))','setup> selecting B minimum number ',CANDIDATES(NRANDOM),
     &                              ' location ',LOCATIONB(CANDIDATES(NRANDOM)),' for the product set'
                  NEWPFMIN(CANDIDATES(NRANDOM))=NEWPFMIN(CANDIDATES(NRANDOM))+1.0D0
               ENDDO
               PFTOTALB=PFTOTALB+NINT(NRWREACTANT*RWPROB(J1))
            ENDDO
            PFTOTALB=LOG(PFTOTALB) ! partition functions are stored as log's
            NCOUNT=0
            DO J1=1,NMINB
               IF (NEWPFMIN(J1).GT.0.0D0) THEN
                  NCOUNT=NCOUNT+1
                  LOCATIONB(NCOUNT)=LOCATIONB(J1)
                  PFMIN(LOCATIONB(NCOUNT))=LOG(NEWPFMIN(J1))
                  PRINT '(A,I8,A,G20.10)','setup> relative weight for reactant minimum ',LOCATIONB(NCOUNT),' is ',
     &                        EXP(PFMIN(LOCATIONB(NCOUNT))-PFTOTALB)
               ENDIF
            ENDDO
            NMINB=NCOUNT
            PRINT '(A,I8,A)','setup> there are now ',NMINB,' minima of type B'
         ELSE
            ALLOCATE(NEWPFMIN(NMINA))
            NEWPFMIN(1:NMINA)=0.0D0
!
!  Select NRWREACTANT minima from the A set according to the required weights in RWPROB
!
            PFTOTALA=0.0D0
            DO J1=1,NRWBINS
               IF (NINT(NRWREACTANT*RWPROB(J1)).EQ.0) CYCLE
               NCOUNT=0
               DO J2=1,NMINA ! identify minima in the required energy range
                  IF ((EMIN(LOCATIONA(J2)).GE.RWEMIN+(J1-1)*RWBINWIDTH).AND.
     &                (EMIN(LOCATIONA(J2)).LE.RWEMIN+J1*RWBINWIDTH)) THEN
                      NCOUNT=NCOUNT+1
                      CANDIDATES(NCOUNT)=J2
                  ENDIF
               ENDDO
               PRINT '(3(A,I8),A,G20.10)','setup> number of A minima in energy bin ',J1,' is ',NCOUNT,' number needed=',
     &                           NINT(NRWREACTANT*RWPROB(J1)),' probability=',RWPROB(J1)
               IF (NCOUNT.EQ.0) STOP
               DO J2=1,NINT(NRWREACTANT*RWPROB(J1))
                  NRANDOM=NCOUNT*DPRAND()+1
                  PRINT '(3(A,I8))','setup> selecting A minimum number ',CANDIDATES(NRANDOM),
     &                              ' location ',LOCATIONA(CANDIDATES(NRANDOM)),' for the product set'
                  NEWPFMIN(CANDIDATES(NRANDOM))=NEWPFMIN(CANDIDATES(NRANDOM))+1.0D0
               ENDDO
               PFTOTALA=PFTOTALA+NINT(NRWREACTANT*RWPROB(J1))
            ENDDO
            PFTOTALA=LOG(PFTOTALA) ! partition functions are stored as log's
            NCOUNT=0
            DO J1=1,NMINA
               IF (NEWPFMIN(J1).GT.0.0D0) THEN
                  NCOUNT=NCOUNT+1
                  LOCATIONA(NCOUNT)=LOCATIONA(J1)
                  PFMIN(LOCATIONA(NCOUNT))=LOG(NEWPFMIN(J1))
                  PRINT '(A,I8,A,G20.10)','setup> relative weight for reactant minimum ',LOCATIONA(NCOUNT),' is ',
     &                        EXP(PFMIN(LOCATIONA(NCOUNT))-PFTOTALA)
               ENDIF
            ENDDO
            NMINA=NCOUNT
            PRINT '(A,I8,A)','setup> there are now ',NMINA,' minima of type A'
         ENDIF
         DEALLOCATE(NEWPFMIN,CANDIDATES)
         ! }}}
      ENDIF


