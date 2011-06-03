
      CALL CPU_TIME(TNEW)

IF (NATTEMPT.GT.0) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent iterating committor probability=',TPFOLD,' s'
      IF (TFOLDT) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent iterating waiting times=',TTFOLD,' s'
      IF (GTT) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent in GT                          =',TGT,' s'
      IF (NGTT) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent in NGT                          =',TGT,' s'
      IF (DIJKSTRAT) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent in Dijkstra                    =',TDIJKSTRA,' s'
      IF (KSHORTESTPATHST) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent in kshortestpaths              =',TKSHORTESTPATHS,' s'
      IF (CONNECTREGIONT) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent in connectdist                 =',TCONNECTDIST,' s'
      WRITE(*,'(A,G15.5,A)')  'main> total CPU time spent in PATHSAMPLE            =',TNEW-TINIT,' s'

