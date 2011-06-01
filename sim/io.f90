

      if (LBFGST) then
         WRITE(LFH,'(A)') 'Nocedal LBFGS minimization'
         WRITE(LFH,'(A,I6)') 'Number of updates before reset in LBFGS=',MUPDATE
         WRITE(LFH,'(A,F20.10)') 'Maximum step size=',MAXBFGS
         WRITE(LFH,'(A,G12.4)') 'Guess for initial diagonal elements in LBFGS=',DGUESS
      ENDIF

      WRITE(LFH,'(A,F15.10)') 'Final quench tolerance for RMS gradient ',CQMAX
      WRITE(LFH,'(A,F15.10)') 'Energy difference criterion for minima=',ECONV
      WRITE(LFH,'(A,I5,A,I5)') 'Maximum number of iterations: sloppy quenches ',MAXIT,' final quenches ',MAXIT2

      IF (DEBUG) THEN
         WRITE(MYUNIT,160) 
160      FORMAT('Debug printing is on')
      ENDIF
       
      WRITE(MYUNIT, '(A,G20.10)') 'Maximum allowed energy rise during a minimisation=',MAXERISE

      IF (TARGET) THEN
         WRITE(MYUNIT,'(A)',ADVANCE='NO') 'Target energies: '
         WRITE(MYUNIT,'(F20.10)',ADVANCE='NO') (TARGETS(J1),J1=1,NTARGETS)
         WRITE(MYUNIT,'(A)') ' '
      ENDIF
