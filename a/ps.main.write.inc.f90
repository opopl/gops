
WRITE(*,'(A)') '*************************************************************************************************'
     ! write{{{
      IF (.NOT.DEBUG) WRITE(*,'(A)') 'debug printing is OFF'
      IF (DEBUG) WRITE(*,'(A)') 'debug printing is ON'
      WRITE(*,'(A,I6,A)') 'Running on ',NCPU,' processors'
      IF (TRIM(ADJUSTL(COPYFILES)).NE.'') 
     &  PRINT '(2A)','The following additional files will be copied to distributed nodes: ',TRIM(ADJUSTL(COPYFILES))
      IF (ADDPT) WRITE(*,'(A)') 'all permutational isomers with different moments of inertia will be enumerated automatically'
      IF (NTAG.GT.0) THEN
         PRINT '(I6,A)',NTAG,' tagged atoms with mass factors:'
         PRINT '(I6,F12.2)',(TAGNUM(J1),TAGFAC(TAGNUM(J1)),J1=1,NTAG)
      ENDIF

      IF (TWOD) WRITE(*,'(A)') 'two-dimensional system'
      IF (BULKT) THEN
         WRITE(*,'(A)') 'bulk system'
         WRITE(*,'(A,3G20.10)') 'x, y, z box lengths: ',BOXLX, BOXLY, BOXLZ
      ENDIF
      IF (.NOT.BULKT) WRITE(*,'(A)') 'finite system (not bulk)'
      WRITE(*,'(A,A2)') 'OPTIM system symbol: ',ZSYM
      WRITE(*,'(A,I5)') 'number of atoms=',NATOMS
      IF (FREEZE) WRITE(*,'(A,I5)') 'number of frozen atoms=',NFREEZE
      IF (AMHT) WRITE(*,'(A,I6)') 'AMH potential, number of glycines=',NGLY
      IF (NTAG.GT.0) WRITE(*,'(A,I5,A,F12.3)') 'tagged atom ',NTAG,' has mass ',TAGMASS
      WRITE(*,'(A,I6)') 'minimum number of connections for minima=',CONNECTIONS
      WRITE(*,'(A,I6)') 'maximum number of ts searches to achieve this minimum value=',MAXTSATTEMPTS
      WRITE(*,'(A,I6)') 'random number seed=',ISEED
      CALL RANDOM_SEED(ISEED)
      WRITE(*,'(A,F15.5)') 'initial perturbation parameter single-ended transition state searches=',PERTVALUE
      WRITE(*,'(A,E20.10)') 'energy difference  criterion for distinguishing stationary points=',EDIFFTOL
      WRITE(*,'(A,E20.10)') 'distance criterion for distinguishing stationary points=',GEOMDIFFTOL
      WRITE(*,'(A,E20.10)') 'inertia difference criterion for setup check=',IDIFFTOL
      IF (TFOLDT) THEN
         PRINT '(A)','Calculation of MFPT and rate constants using Gauss-Seidel/successive overrelaxation iteration'
         PRINT '(A,F15.1)','Maximum steps=',NTFOLD
         PRINT '(A,G20.10)','Convergence condition=',TFOLDTHRESH
         PRINT '(A,G20.10)','Relaxation factor (omega=1 is Gauss-Seidel)=',TOMEGA
      ENDIF
      IF ((.NOT.GTT).AND.(.NOT.GT2T).AND.(.NOT.NGTT)) THEN
         IF (DIRECTION.EQ.'AB') WRITE(*,'(A)') 'sampling direction is a<-b'
         IF (DIRECTION.EQ.'BA') WRITE(*,'(A)') 'sampling direction is b<-a'
         WRITE(*,'(A,A80)') 'OPTIM executable: ',EXEC
      ENDIF
      IF (REMOVEUNCONNECTEDT) THEN
         WRITE(*,'(A,I5,A)') 'Rewriting database removing minima with ',NCONNMIN,' connections or fewer'
         WRITE(*,'(A)') 'and stationary points unconnected to ' // TRIM(ADJUSTL(UNCONNECTEDS))
      ENDIF
      IF (ENSEMBLE.EQ.'T') WRITE(*,'(A,F15.5)') 'Temperature=',TEMPERATURE
      IF (ENSEMBLE.EQ.'E') WRITE(*,'(A,F15.5)') 'Total energy=',TOTALE
      IF (NOFRQS) WRITE(*,'(A)') 'No frequency information should be present in info files'
      IF (KMCT) THEN
         WRITE(*,'(A,F20.1,A,F12.4)') 'Running KMC on database. Number of KMC runs=',NKMCCYCLES,
     1                             ' pair renorm threshold=',PAIRTHRESH
      ELSE IF (DIJKSTRAT.AND.(.NOT.DIJPAIRT)) THEN
         WRITE(*,'(A)') 'Performing Dijkstra analysis to find the most important paths'
         IF (KSHORTESTPATHST) THEN
             WRITE(*,'(A)') 'ERROR: KSHORTESTPATHS and DIJKSTRA are incompatible; stopping.'
             STOP
         ENDIF
      ELSE IF (KSHORTESTPATHST) THEN
         WRITE(*,'(A)') 'Performing k-th shortest path analysis to find the most important paths'
         IF (DIJKSTRAT) THEN
             WRITE(*,'(A)') 'ERROR: KSHORTESTPATHS and DIJKSTRA are incompatible; stopping.'
             STOP
         ENDIF
      ELSE IF (KMCCOMMITT) THEN
         WRITE(*,'(A)')        'Using committor probability/KMC routine:'
         WRITE(*,'(A,F20.1)')     '      KMC cycles per starting minimum=                 ',NKMCCYCLES
         WRITE(*,'(A,F20.1)')  '      Maximum number of steps in a KMC run=            ',MAXBREAK
         WRITE(*,'(A,G20.10)') '      Convergence condition on committor probabilities=',PABCONV
         WRITE(*,'(A,G20.10)') '      Successive overrelaxation parameter omega=       ',OMEGA
         WRITE(*,'(A,G20.10)') '      Threhold for leapfrog KMC moves=                 ',PAIRTHRESH
      ELSE IF (NGTT) THEN
          WRITE(*,'(A,I6,A,F12.4)') 'Running new graph theory rate calculation on database'
      ELSE IF (ARNOLDIT) THEN
          WRITE(*,'(A)') 'Eigenvalues of the transition matrix from weighted Arnoldi subspace method'
      ELSE IF (DIAGT) THEN
          WRITE(*,'(A,I6,A)') 'Explicit diagonalisation of the transition matrix to find ',NDIAGEIG,' eigenvalues'
          WRITE(*,'(A,G12.4)') 'All transition matrix elements will be scaled by ',DIAGSCALE
      ELSE IF (GTT.or.GT2T) THEN
          WRITE(*,'(A,I6,A,F12.4)') 'Running graph theory rate calculation on database'
      ELSE 
         WRITE(*,'(A,I8)') 'Number of cycles for connections between local minima=',NATTEMPT
         IF (PFOLDINT.NE.0) THEN
            PRINT '(A,I6)','Interval between calls to Pfold routine=',PFOLDINT 
            PRINT '(A,I6,A,G20.10)','SOR iterations per call=',NPFOLD,' with omega=',OMEGA 
         ELSE
            PRINT '(A,I6)','No Pfold iterations'
         ENDIF
      ENDIF
      IF (DIJKSTRAT.AND.(MAXDOWNBARRIER.LT.1.0D99)) THEN
         PRINT '(A,G20.10)','Downhill barrier check in Dijkstra, threshold=',MAXDOWNBARRIER
      ENDIF
      
      IF (FRICTIONT) WRITE(*,'(A,G20.10)') 'Phenomenological friction coefficient in frequency units=',GAMMAFRICTION
      IF (REGROUPT) WRITE(*,'(A,F15.5)') 'Reclassifying A/B/I minima for threshold energy ',REGROUPTHRESH
      IF (REGROUPRATET) THEN
         WRITE(*,'(A,G20.10)') 'Grouping and calculating free energies for threshold rate ',REGROUPRATETHRESH
         WRITE(*,'(A,G20.10)') 'Plank constant set to ',PLANCK
         IF (CHARMMT.AND.(PLANCK.EQ.1.0D0)) THEN
            PRINT '(2(A,G20.10))','ERROR - regrouping with CHARMM and Planck constant=',PLANCK,' resetting to',9.546D-14
            PLANCK=9.546D-14
         ENDIF
      ENDIF
      IF (REGROUPPET) THEN
         WRITE(*,'(A,G20.10)') 'Regrouping A/B/I minima for superbasin PE ',REGROUPPETHRESH
      ENDIF
      IF (REGROUPFREET) THEN
         WRITE(*,'(A,G20.10)') 'Grouping and calculating free energies for threshold barrier ',REGROUPFREETHRESH
         WRITE(*,'(A,G20.10)') 'Plank constant set to ',PLANCK
         IF (CHARMMT.AND.(PLANCK.EQ.1.0D0)) THEN
            PRINT '(2(A,G20.10))','ERROR - regrouping with CHARMM and Planck constant=',PLANCK,' resetting to',9.546D-14
            PLANCK=9.546D-14
         ENDIF
      ENDIF
      IF (REGROUPFREEABT) THEN
         WRITE(*,'(A,G20.10)') 'Grouping and calculating free energies for threshold barrier ',REGROUPFREETHRESH
         WRITE(*,'(A,G20.10)') 'Free energy groups will not be used except to expand the A and B sets'
         WRITE(*,'(A,G20.10)') 'Plank constant set to ',PLANCK
         IF (CHARMMT.AND.(PLANCK.EQ.1.0D0)) THEN
            PRINT '(A,G20.10)','ERROR - regrouping with CHARMM and Planck constant=',PLANCK
            PRINT '(2(A,G20.10))','ERROR - regrouping with CHARMM and Planck constant=',PLANCK,' resetting to',9.546D-14
         ENDIF
      ENDIF
      WRITE(*,'(A,I5,A)') 'Rate calculation will iteratively remove minima with ',NCONNMIN,' connections or fewer'
      IF (STARTFROMPATH) THEN
         IF (STARTTRIPLES) THEN
            WRITE(*,'(A)') 'Reading initial pathway from ' // TRIM(ADJUSTL(PATHNAME)) //' in min-sad-min format'
         ELSE 
            WRITE(*,'(A)') 'Reading initial pathway from ' // TRIM(ADJUSTL(PATHNAME)) //' in old format'
         ENDIF
         PRINT '(2(A,I8))','The A minimum is number ',STARTMINA,' and the B minimum is ',STARTMINB
      ENDIF
      IF (ADDPATH) THEN
         IF (ADDTRIPLES) THEN
            WRITE(*,'(A)') 'Adding stationary points from ' // TRIM(ADJUSTL(PATHNAME)) //' in min-sad-min format'
         ELSE
            WRITE(*,'(A)') 'Adding stationary points from ' // TRIM(ADJUSTL(PATHNAME)) // ' in old format'
         ENDIF
      ENDIF
      IF (TRIPLES) WRITE(*,'(A)') 'path.info files will be read as min-sad-min triples'
      IF (CONNECTREGIONT) THEN
         PRINT '(A,G12.3,A,2I8)','Pairs of minima will be chosen based upon pairs separated by less than ',  
     &               CONNECTDIST,' starting from minima ',CONNECTMIN1,CONNECTMIN2
      ELSEIF (DIJINITT) THEN
         WRITE(*,'(A)') 'Performing Dijkstra analysis to perform initial connection'
         IF (INDEXCOSTFUNCTION) THEN
            WRITE(*,'(A)') 'Using index of minima for edge weights'
         ELSEIF (EXPCOSTFUNCTION) THEN
            WRITE(*,'(A)') 'Using exponential distance for edge weights'
         ELSEIF (.NOT.INTCONSTRAINTT) THEN
            WRITE(*,'(A,I4,A)') 'Using distance to the power ',CostFunctionPower,' for edge weights'
         ENDIF
         WRITE(*,'(A,G20.10)') 'Transition states will be rejected above threshold energy ',TSTHRESH
      ELSEIF (DIJINITFLYT) THEN
         WRITE(*,'(A)') 'Performing Dijkstra analysis to perform initial connection with weights calculated on-the-fly'
         IF (EXPCOSTFUNCTION) THEN
            WRITE(*,'(A)') 'Using exponential distance for edge weights'
         ELSEIF (.NOT.INTCONSTRAINTT) THEN
            WRITE(*,'(A,I4,A)') 'Using distance to the power ',CostFunctionPower,' for edge weights'
         ENDIF
         WRITE(*,'(A,G20.10)') 'Transition states will be rejected above threshold energy ',TSTHRESH
      ELSEIF (DIJPAIRT) THEN
         WRITE(*,'(A)') 'Performing Dijkstra analysis to generate connnection pairs'
         IF (BARRIERSORT) THEN
            WRITE(*,'(A)') 'Pairs will be sorted according to barrier height'
         ELSE
            WRITE(*,'(A)') 'Pairs will be sorted according to energy of the minimum'
         ENDIF
!     ELSEIF (PTAUT) THEN
!        PRINT '(A,I8,A)','Minima will be chosen for connection to minima on the fastest path based on p^eq * waiting time'
      ELSEIF (FREEPAIRT) THEN
         PRINT '(A,I8,A)','Pairs of minima will be chosen based on free energy barriers between regrouped minima'
         PRINT '(3(A,G20.10))','Regrouping barrier threshold=',REGROUPFREETHRESH,' with minima above ',FREETHRESH,' excluded'
         IF (NPAIRFRQ.LT.1) THEN
            PRINT '(A,F12.3)','Pairs will not be recalculated unless we run out of possibilities'
         ELSE
            PRINT '(A,I6,A)','Pairs will be recalculated every ',NPAIRFRQ,' cycles'
         ENDIF
         IF (CHARMMT.AND.(PLANCK.EQ.1.0D0)) THEN
            PRINT '(A,G20.10)','ERROR - regrouping with CHARMM and Planck constant=',PLANCK
            STOP
         ENDIF
      ELSEIF (SHORTCUTT) THEN
         IF (BARRIERSHORT) THEN
            PRINT '(A,I8,A)','Pairs of minima will be chosen for minima separated by fewer than',MINSEP,  
     &                  ' steps on either side of the largest barriers in the best SS path'
         ELSEIF (RATESHORT) THEN
            PRINT '(A,I8,A)','Pairs of minima will be chosen for minima separated by fewer than',MINSEP,  
     &                  ' steps on either side of the smallest rate constants in the best SS path'
         ELSE
            PRINT '(A,I8,A)','Pairs of minima will be chosen for the closest minima separated by at least ',MINSEP,  
     &                  ' steps in the best SS path identified by Dijkstra'
         ENDIF
         IF (NPAIRFRQ.LT.1) THEN
            PRINT '(A,F12.3)','Pairs will not be recalculated unless we run out of possibilities'
         ELSE
            PRINT '(A,I6,A)','Pairs will be recalculated every ',NPAIRFRQ,' cycles'
         ENDIF
      ELSEIF (UNTRAPT) THEN
         PRINT '(A,I8,A)','Pairs of minima will be chosen according to their lowest barrier to a minimum in a product superbasin'
         PRINT '(A,F15.5,A)','Connections will not be attempted to minima higher than ',EUNTRAPTHRESH, 
     &                      ' energy units above the first minimum'
         IF (NPAIRFRQ.LT.1) THEN
            PRINT '(A,F12.3)','Pairs will not be recalculated unless we run out of possibilities'
         ELSE
            PRINT '(A,I6,A)','Pairs will be recalculated every ',NPAIRFRQ,' cycles'
         ENDIF
      ELSEIF (NATTEMPT.GT.0) THEN
         PRINT '(A,F12.3)','Pairs of minima will be chosen based upon a minimum committor probability difference of ',PSCALE
         IF (NPAIRFRQ.LT.1) THEN
            PRINT '(A,F12.3)','Pairs will not be recalculated unless we run out of possibilities'
         ELSE
            PRINT '(A,I6,A)','Pairs will be recalculated every ',NPAIRFRQ,' cycles'
         ENDIF
      ELSEIF (USEPAIRST) THEN
         PRINT '(A,A)','Pairs of minima will be chosen from the sequence in file ',TRIM(ADJUSTL(USEPAIRSFILE))
      ENDIF
      IF (DUMMYTST) THEN
         IF (BHINTERPT) PRINT '(A)','Dummy ts entries will be created for closest pairs of minima found via basin-hopping.'
         IF (BISECTT) PRINT '(A)','Dummy ts entries will be created for pairs of minima found by bisection.'
      ENDIF
      IF (CLOSEFILEST) PRINT '(A)','The min.data and ts.data files will not be kept open'
      IF (INTLJT) THEN
         PRINT '(A,F15.5)',   'Using interpLJ potential metric'
         PRINT '(A,F20.10)',  'minimum distance difference for internal minimum=',INTLJDEL
         PRINT '(A,F20.10)',  'multiplying factor for internal minimum penalty function=',INTLJEPS
      ENDIF
      IF (INTCONSTRAINTT) THEN
         PRINT '(A,F15.5)',   'Using constraint potential metric with absolute distance change tolerance ',INTCONSTRAINTTOL
         PRINT '(A,F15.5)',   '   constraint spring constant=',INTCONSTRAINTDEL
         PRINT '(2(A,F15.5))','   repulsion factor between unconstrained atoms=',INTCONSTRAINTREP
         PRINT '(A,F15.5)',   '   cutoff for repulsion will be the minimum of ',INTCONSTRAINREPCUT
         PRINT '(A)',         '   and the shortest distance in the end points'
         PRINT '(A,I6)',      '   maximum separation of atoms in sequence for constraint=',INTCONSEP
         PRINT '(A,I6)',      '   minimum separation of atoms in sequence for repulsion=',INTREPSEP
         PRINT '(A,I6)',      '   maximum number of constraints per atom=',MAXCONUSE
         IF (CHECKCONINT) THEN
            PRINT '(A,I6)',   '   adding terms for constraint internal minima'
         ELSE
            PRINT '(A,I6)',   '   not adding terms for constraint internal minima'
         ENDIF
      ENDIF 
      ! }}}
      WRITE(*,'(A)') '*************************************************************************************************'


