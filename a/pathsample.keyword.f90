!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales                      {{{
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
!                                                                }}}
!op226 <begin>
      SUBROUTINE KEYWORD
      ! Declarations 
      ! {{{
      ! Modules {{{
      USE PORFUNCS
      USE NODES, ONLY: JPN, GETNODES, NNODES
      USE KEY
      USE COMMON
      USE RIGIDBODYMOD, ONLY: CAPSOMER, NUMRBTYPES, CAPSOMERDEFS
        ! }}}
      IMPLICIT NONE

      ! Variables {{{
      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, ISTAT, NDUMMY2, LAST
      LOGICAL CAT
      COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR, NERROR, IR, ECHO, LAST, CAT
      INTEGER J1, NDUMMY, J2, J3
      DOUBLE PRECISION DUMMY, DBSIGBB
      LOGICAL END, SKIPBL, CLEAR, ECHO, PERMFILE, RBSYMTEST
      CHARACTER(LEN=20) WW
      CHARACTER(LEN=80) RWENERGYFILE
      CHARACTER(LEN=20) WORD
      CHARACTER(LEN=1) LOYNO
      CHARACTER(LEN=9) UNSTRING

      TYPE (CAPSOMER), ALLOCATABLE :: TEMPCAPSOMERDEFS(:)
      INTEGER :: MAXRBTYPES
      ! }}}
!}}}
!op226 </begin>
!op226 Initializations {{{ 
!op226 <init>
      DBPT     = .FALSE.
      DBPTDT   = .FALSE.
      MSSTOCKT = .FALSE.
      PAHAT    = .FALSE.
      PATCHYDT = .FALSE.
      STOCKAAT = .FALSE.
      EFIELDT  = .FALSE.
      DEBUG=.FALSE.
      EXTRACTMINT=.FALSE.
      EXTRACTTST=.FALSE.
      ADDPT=.FALSE.
      DUMPGROUPST=.FALSE.
      DSCALE=3.0D0
      PSCALE=0.5D0
      CONNECTIONS=0
      MAXTSATTEMPTS=10
      ISEED=1
      NATTEMPT=0
      PERTVALUE=0.9D0
      PERTMAX=2.0D0*PERTVALUE
      PERTMIN=0.5D0*PERTVALUE
      EDIFFTOL=1.0D-8
      GEOMDIFFTOL=1.0D-1
      IDIFFTOL=1.0D-3
      DIRECTION='AB'
      ENSEMBLE='T'
      TWOD=.FALSE.
      BULKT=.FALSE.
      ANGLEAXIS=.FALSE.
      NUMRBTYPES = 0
      MAXRBTYPES = 0

      AMBERT=.FALSE.
      AMHT=.FALSE.
      AMHQT=.FALSE.
      AMHQCONTT=.FALSE.
      AMHRMSDT=.FALSE.
      AMHRELQT=.FALSE.
      AMH_RELCOT=.FALSE.
      AMHALLATOMMINT=.FALSE.
      AMHALLATOMTST=.FALSE.
      NOFRQS=.FALSE.

      JPN=1
      NNODES=1
      NATOMS=-1
      CHARMMT=.FALSE.
      STARTFROMPATH=.FALSE.
      READMINT=.FALSE.
      ADDPATH=.FALSE.
      KMCCOMMITT=.FALSE.
      GTT=.FALSE.
      NGTT=.FALSE.
      NGTDISCONNECTALL=.FALSE.
      NGTSWITCH=0.3D0
      NGTSIZE=11000
      NGTCRSWITCH=2.0D0
C     Sem: begin GT2 controls
      GT2T=.FALSE.
      GT2SPARSE=.TRUE.
      GT2SWITCH=.TRUE.
      GT2DisconnectSources=.TRUE.
      GT2ALTPBB=.TRUE.
      GT2RESCALE=.FALSE.
      GT2NORMALISE=.FALSE.
      GT2RSWITCH=0.08D0
      GT2PTOL=1.0D-5
C     Sem: end GT2 controls
      GTINT=0
      REGROUPT=.FALSE.
      REGROUPTHRESH=-1.0D100
      REGROUPRATET=.FALSE.
      REGROUPRATETHRESH=1.0D100
      REGROUPPET=.FALSE.
      REGROUPPETHRESH=-1.0D100
      REGROUPFREET=.FALSE.
      REGROUPFREEABT=.FALSE.
      REGROUPFREETHRESH=-1.0D100
      PABCONV=1.0D-8
      OMEGA=1.0D0 ! GAUSS-SEIDEL ITERATION
      KMCT=.FALSE.
      NCONNMIN=0
      MAXBREAK=1.0D6
      NKMCCYCLES=100.0D0
      PAIRTHRESH=1.0D0
      UNRST=.FALSE.
      NINTS=-1
      NCPU=1
      UNSTRING='UNDEFINED'
      ZSYM='UNKNOWN'
      EXEC='UNDEFINED'
      EXECGMIN='UNDEFINED'
      PATHNAME='UNDEFINED'
      NOPOINTS=.FALSE.
      TRIPLES=.FALSE.
      STARTTRIPLES=.FALSE.
      ADDTRIPLES=.FALSE.
      DIJKSTRAT=.FALSE.
      DIJKSTRAWAITT=.FALSE.
      DIJPAIRT=.FALSE.
      BARRIERSORT=.FALSE.
      BARRIERSHORT=.FALSE.
      RATESHORT=.FALSE.
      DIJINITT=.FALSE.
      DIJINITFLYT=.FALSE.
      DIJINITSTARTT=.FALSE.
      DIJINITCONTT=.FALSE.
      TSTHRESH=HUGE(1.0D0)
      MAXBARRIER=HUGE(1.0D0)
      MAXDOWNBARRIER=HUGE(1.0D0)
      COSTFUNCTIONPOWER=1
      EXPCOSTFUNCTION=.FALSE.
      INDEXCOSTFUNCTION=.FALSE.
      COPYFILES=''
      COPYOPTIMT=.FALSE.
      NPFOLD=0
      TFOLDT=.FALSE.
      NTFOLD=1.0D4 ! real not integer !
      TOMEGA=1.0D0
      TFOLDTHRESH=1.0D-3
      MAXATTEMPT=1
      CALCORDERT=.FALSE.
      CONNECTREGIONT=.FALSE.
      CONNECTMIN1=1
      CONNECTMIN2=1
      CONNECTDIST=1.0D10
      SHORTCUTT=.FALSE.
!     PTAUT=.FALSE.
      NPAIRFRQ=0
      MERGEDBT=.FALSE.
      UNTRAPT=.FALSE.
      EUNTRAPTHRESH=1.0D100
      FREEPAIRT=.FALSE.
      PERMDIST=.FALSE.
      PERMISOMER=.FALSE.
      NTAG=0
      TAGT=.FALSE.
      FREEZE=.FALSE.
      NFREEZE=0
      PLANCK=1.0D0
      DUMMYRUNT=.FALSE.
      DUMMYTST=.FALSE.
      REWEIGHTT=.FALSE.
      KSHORTESTPATHST = .FALSE.
      KSHORT_FULL_PRINTT = .FALSE.
      NPATHS = 0
      BHINTERPT=.FALSE.
      BHACCREJ=0.5D0
      BHSTEPSIZE=0.4D0
      BHCONV=0.01D0
      BHSTEPS=1000
      BHTEMP=1.0D0
      BHK=1.0D0
      ICINTERPT=.FALSE.

      USEPAIRST=.FALSE.
      LOWESTFRQT=.FALSE.
      IMFRQT=.FALSE.
      EVCUT=2.0D-6

      BISECTT=.FALSE.
      DIAGT=.FALSE.
      ARNOLDIT=.FALSE.
      SLURMT=.FALSE.
      CVT=.FALSE.
      DOST=.FALSE.
      CHECKCONNECTIONST=.FALSE.
      CLOSEFILEST=.FALSE.
      PULLT=.FALSE.
!     DC430 >
      RBAAT  = .FALSE.
      RBSYMT = .FALSE.
      FRICTIONT=.FALSE.
      GAMMAFRICTION=0.0D0
      REMOVEUNCONNECTEDT=.FALSE.
      UNCONNECTEDS='AB'
!
! Constraint potential
!
      INTERPCOSTFUNCTION=.FALSE.
      INTCONSTRAINTT=.FALSE.
      INTCONSTRAINTTOL=0.1D0
      INTCONSTRAINTDEL=1.0D5
      INTCONSTRAINTREP=1.0D0
      INTCONSTRAINREPCUT=20.0D0
      INTREPSEP=0
      INTCONSEP=10000
      CHECKCONINT=.FALSE.
      MAXCONUSE=3
!
! LJ interpolation potential parameters
!
      INTLJT=.FALSE.
      INTLJDEL=0.1D0
      INTLJEPS=1.0D0
C
C
C sf344> docking stuff
C
      DOCKT=.FALSE.
      DSTAGE(:)=.TRUE.
!
! SIS epidemiological model
! 
      SIST=.FALSE.
      SMAX=0
      IMAX=0
      POPSS=0
      SISMU=0.0D0
      SISKAPPA=0.0D0
      SISBETA=0.0D0

!op226 </init>
!op226 }}}
!op226 Open and read pathdata {{{
!op226 <read>
C
C SAT: Now, there it is! Unit 5 is a standard input unit, and we must not use it!
C 5--->959
C
      OPEN (959,FILE='pathdata',STATUS='OLD')
C
C  Read from unit 959
C
      IR=959
190   CALL INPUT(END)
      IF (.NOT. END) THEN
        CALL READU(WORD)
      ENDIF

      IF (END .OR. WORD .EQ. 'STOP') THEN
        CLOSE(IR)
        RETURN
      ENDIF

      IF (WORD.EQ.'    '.OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'
     &                          .OR. WORD .EQ. '\\') THEN 
         GOTO 190

! </read> 
!op226 }}}

!op226 keyword loop {{{
! <kwd>
C
C  Add the minima and transition= states in path.info.<PATHNAME>
C  to an existing database. The end points are NOT assumed to
C  belong to the A and B sets.
C

include pathsample.keyword.addpath.inc.f90
include pathsample.keyword.addperm.inc.f90
include pathsample.keyword.addtriples.inc.f90
include pathsample.keyword.angleaxis.inc.f90
include pathsample.keyword.amber9.inc.f90
include pathsample.keyword.nab.inc.f90
include pathsample.keyword.amh.inc.f90
include pathsample.keyword.amhq.inc.f90
include pathsample.keyword.amhqcont.inc.f90
include pathsample.keyword.amhqrmsd.inc.f90
include pathsample.keyword.amhqrelq.inc.f90
include pathsample.keyword.amhq_relco.inc.f90
            

C
C  Output AMH all-atom MIN .
C
      ELSE IF (WORD.EQ.'AMHALLATOMMIN') THEN
         AMHALLATOMMINT=.TRUE.
         PRINT '(A,I6,A,I6)','keywords> output AMH all atom MIN structures'
C
C  Output AMH all-atom TS .
C
      ELSE IF (WORD.EQ.'AMHALLATOMTS') THEN
         AMHALLATOMTST=.TRUE.
         PRINT '(A,I6,A,I6)','keywords> output AMH all atom TS structures'

C
C  Rates from eigenvalues obtained with weighted Arnoldi subspace method.
C
      ELSE IF (WORD.EQ.'ARNOLDI') THEN
         ARNOLDIT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
C
C  Parameters for OPTIM basin-hopping interpolation jobs.
C
      ELSE IF (WORD.EQ.'BHINTERP') THEN
         BHINTERPT=.TRUE.
         CALL READF(BHDISTTHRESH)
         CALL READF(BHMAXENERGY)
         CALL READI(BHSTEPS)
         CALL READF(BHCONV)
         CALL READF(BHTEMP)
         CALL READF(BHSTEPSIZE)
         CALL READF(BHACCREJ)
         CALL READF(BHK)
         CALL READF(BHSFRAC)
         IF (NITEMS.GT.10) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ICINTERP') ICINTERPT=.TRUE.
C
C  Parameters for OPTIM bisection interpolation jobs.
C
      ELSE IF (WORD.EQ.'BISECT') THEN
         BISECTT=.TRUE.
         CALL READF(BISECTMINDIST)
         CALL READF(BISECTMAXENERGY)
         CALL READI(BISECTSTEPS)
         CALL READI(BISECTMAXATTEMPTS)
         IF (NITEMS.GT.5) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ICINTERP') ICINTERPT=.TRUE.
C
C  Bulk with periodic boundary conditions.
C
      ELSE IF (WORD.EQ.'BULK') THEN
         BULKT=.TRUE.
         CALL READF(BOXLX)
         CALL READF(BOXLY)
         CALL READF(BOXLZ)
C
C  Calculate order parameter for all stationary points and exit.
C  The order parameter routine will be system specific and must
C  replace the dummy routine provided.
C  Minima and transition states will be listed in separate files
C  for points with order parameters above/below the specified threshold,
C  with default value 0.0.
C
      ELSE IF (WORD.EQ.'CALCORDER') THEN
         CALCORDERT=.TRUE.
C
C Capsomere definition
C
      ELSE IF (WORD.EQ.'CAPSID') THEN
         ANGLEAXIS=.TRUE.
         RIGIDBODY=.TRUE.
         IF (.NOT.ALLOCATED(CAPSOMERDEFS)) THEN
            MAXRBTYPES = 10
            ALLOCATE(CAPSOMERDEFS(MAXRBTYPES))

         ELSE IF (NUMRBTYPES.EQ.MAXRBTYPES) THEN
            ALLOCATE(TEMPCAPSOMERDEFS(MAXRBTYPES))
            TEMPCAPSOMERDEFS = CAPSOMERDEFS

            MAXRBTYPES = MAXRBTYPES + 10
            DEALLOCATE(CAPSOMERDEFS)
            ALLOCATE(CAPSOMERDEFS(MAXRBTYPES))

            CAPSOMERDEFS(1:NUMRBTYPES) = TEMPCAPSOMERDEFS
            DEALLOCATE(TEMPCAPSOMERDEFS)
         ENDIF

         NUMRBTYPES = NUMRBTYPES + 1

         CAPSOMERDEFS(NUMRBTYPES)%NBASALSITES = 5

         CALL READF(CAPSOMERDEFS(NUMRBTYPES)%RHO)
         CALL READF(CAPSOMERDEFS(NUMRBTYPES)%EPSILONREP)
         CALL READF(CAPSOMERDEFS(NUMRBTYPES)%RADIUS)
         IF (NITEMS.GT.4) THEN
            CALL READF(CAPSOMERDEFS(NUMRBTYPES)%HEIGHT)
            IF (NITEMS.GT.5) CALL READI(CAPSOMERDEFS(NUMRBTYPES)%NBASALSITES)
         ELSE
            CAPSOMERDEFS(NUMRBTYPES)%HEIGHT = 0.5D0*CAPSOMERDEFS(NUMRBTYPES)%RADIUS
         ENDIF
C
C  Set CHARMM potential. NDIHE is the number of dihedral angles in the system, which
C is used in perturb and tssearch to perturb the system randomly
C
      ELSE IF (WORD.EQ.'CHARMM') THEN
         CHARMMT=.TRUE.
         CALL READI(NDIHE)
C
C  Check minimum number of connections for each minimum at startup.
C
      ELSE IF (WORD.EQ.'CHECKCONNECTIONS') THEN
         CHECKCONNECTIONST=.TRUE.
C
C  Check for internal minimum in constraint terms for INTCONSTRAINT
C
      ELSE IF (WORD.EQ.'CONINT') THEN
         CHECKCONINT=.TRUE.
C
C  Close and open all files as needed. The objective is to work around
C  the misinteraction of nfs with the Linux kernel, which causes random
C  control characters to be written to open nfs mounted files.
C
      ELSE IF (WORD.EQ.'CLOSEFILES') THEN
         CLOSEFILEST=.TRUE.
C
C  Minimum number of connections for each minimum.
C
      ELSE IF (WORD.EQ.'CONNECTIONS') THEN
         CALL READI(CONNECTIONS)
         IF (NITEMS.GT.2) CALL READI(MAXTSATTEMPTS)
C
C  Specify that the database is grown by attempting all connections between
C  known minima within a distance cutoff.
C
      ELSE IF (WORD.EQ.'CONNECTREGION') THEN
         CONNECTREGIONT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(CONNECTMIN1)
         IF (NITEMS.GT.2) CALL READI(CONNECTMIN2)
         IF (NITEMS.GT.3) CALL READF(CONNECTDIST)
C
C  Specify files to be copied to distributed nodes in addition to the default
C  odata.<pid> and finish.<pid> files.
C
      ELSE IF (WORD.EQ.'COPYFILES') THEN
C        PRINT '(A,I8)','NITEMS=',NITEMS
         DO J1=2,NITEMS
            CALL READA(WW)
            WRITE(COPYFILES,'(3A)') TRIM(ADJUSTL(COPYFILES)),' ',TRIM(ADJUSTL(WW))
C           PRINT '(3A)','WW,COPYFILES: ',WW,COPYFILES
         ENDDO
C
C  All OPTIM output files will be copied from distributed nodes and not deleted. 
C  This also occurs if DEBUG is .TRUE.
C
      ELSE IF (WORD.EQ.'COPYOPTIM') THEN
         COPYOPTIMT=.TRUE.
C
C  Use NCPU's by starting up to NCPU's OPTIM jobs.
C
      ELSE IF (WORD.EQ.'CPUS') THEN
         CALL READI(NCPU)
      ELSE IF (WORD.EQ.'CV') THEN
         CVT=.TRUE.
         CVTMIN=1.0D0
         CVTMAX=2.0D0
         CVTINC=0.1D0
         CALL READF(CVTMIN)
         CALL READF(CVTMAX)
         CALL READF(CVTINC)
C
C  "JOBSPERNODE JPN" specifies parallel run with JPN jobs per node. Mutually exclusive with CPUS keyword.
C  This requests to read node file. Number of nodes and nodefile must be
C  specified on the command line for this to work. If no command line args are
C  given "JOBSPERNODE JPN" acts same as "CPUS NCPUS".
C
C
C  Maximum number of connect cycles in CYCLE.
C
      ELSE IF (WORD.EQ.'CYCLES') THEN
         CALL READI(NATTEMPT)

      ELSE IF (WORD.EQ.'DB') THEN
! Set the EFIELDT separately, if required
         DBPT = .TRUE.
         CALL READF(DBSIGBB)
         NRBSITES = 3
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL DEFDUM(DBSIGBB)

      ELSE IF (WORD.EQ.'DBTD') THEN
! Set the EFIELDT separately, if required
         DBPTDT = .TRUE.
         CALL READF(DBSIGBB)
         NRBSITES = 3
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL DEFDUM(DBSIGBB)
C
C  Turn on debug printing.
C
      ELSE IF (WORD.EQ.'DEBUG') THEN
         DEBUG=.TRUE.
C
C  Rates from direct diagonalisation.
C
      ELSE IF (WORD.EQ.'DIAG') THEN
         DIAGT=.TRUE.
         NDIAGEIG=10
         DIAGSCALE=1.0D0
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
         IF (NITEMS.GT.2) CALL READI(NDIAGEIG)
         IF (NITEMS.GT.3) CALL READF(DIAGSCALE)
C
C  DIJINIT specifies a Dijkstra analysis to try and construct an initial path.
C  We only need files start and finish containing the end point coordinates for
C  DIJINITSTART - this keyword creates new min.A, min.B, points.min, points.ts,
C  min.data and ts.data files, so it is potentially dangerous!
C  DIJINITCONT specifies an initial path run where the above files already exist.
C  In a DIJINITSTART run C  file min.A is created exist with entries 1 1 and file min.B
C  with entries 1 2. min.data is set up with the two entries for
C  these minima. The coordinates are obtained via odata.start and odata.finish and put in records
C  1 and 2 in the points.min file. 
C
      ELSE IF (WORD.EQ.'DIJINITCONT') THEN
         IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='INDEX') then
                    INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='EXP') then
                    EXPCOSTFUNCTION = .TRUE.
               ELSE IF (WW(1:1) /= ' ') then
                    READ(WW,'(I20)') COSTFUNCTIONPOWER
               ENDIF
         ENDIF
         DIJINITCONTT=.TRUE.
         DIJINITT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
      ELSE IF (WORD.EQ.'DIJINITSTART') THEN
         IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='INDEX') then
                    INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='EXP') then
                    EXPCOSTFUNCTION = .TRUE.
               ELSE IF (WW(1:1) /= ' ') then
                    READ(WW,'(i20)') CostFunctionPower
               ENDIF
         ENDIF
         DIJINITSTARTT=.TRUE.
         DIJINITT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
      ELSE IF (WORD.EQ.'DIJINITCONTFLY') THEN
         IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='INDEX') then
                    INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='EXP') then
                    EXPCOSTFUNCTION = .TRUE.
               ELSE IF (WW(1:1) /= ' ') then
                    READ(WW,'(I20)') COSTFUNCTIONPOWER
               ENDIF
         ENDIF
         DIJINITCONTT=.TRUE.
         DIJINITFLYT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
      ELSE IF (WORD.EQ.'DIJINITSTARTFLY') THEN
         IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='INDEX') then
                    INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='EXP') then
                    EXPCOSTFUNCTION = .TRUE.
               ELSE IF (WW(1:1) /= ' ') then
                    READ(WW,'(i20)') CostFunctionPower
               ENDIF
         ENDIF
         DIJINITSTARTT=.TRUE.
         DIJINITFLYT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
C
C  Dijkstra analysis to find the largest contribution to the SS rate constant
C  from each B (or A) minimum.
C
      ELSE IF (WORD.EQ.'DIJKSTRA') THEN
         DIJKSTRAT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
C
C  Dijkstra analysis to find the path with the lowest sum of waiting times
C  times initial conditional probability.
C
      ELSE IF (WORD.EQ.'DIJKSTRAWAIT') THEN
         DIJKSTRAWAITT=.TRUE.
         DIJKSTRAT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
C
C  DIJPAIR specifies a Dijkstra analysis to propose the pairs of minima to
C  try and connect during the database construction.
C
      ELSE IF (WORD.EQ.'DIJPAIR') THEN
         CALL READI(MAXATTEMPT)
         IF (NITEMS.GT.2) CALL READI(NCONNMIN)
         IF (NITEMS.GT.3) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'BARRIER') BARRIERSORT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
C
C  Direction of sampling.
C
      ELSE IF (WORD.EQ.'DIRECTION') THEN
         CALL READA(DIRECTION)
C
C  Docking routine to calculate the binding free energy of a ligand to a protein
C
      ELSE IF (WORD.EQ.'DOCK') THEN
         DOCKT=.TRUE.
         CALL READI(PARALLEL)
         IF (NITEMS.GT.2) THEN
          DO J1=1,5
            CALL READI(NDUMMY)
            IF(NDUMMY==0) DSTAGE(J1)=.FALSE.
          END DO
         END IF
      ELSE IF (WORD.EQ.'DOS') THEN
         DOST=.TRUE.
         DOSEMIN=1.0D0
         DOSEMAX=2.0D0
         DOSEINC=0.1D0
         CALL READF(DOSEMIN)
         CALL READF(DOSEMAX)
         CALL READF(DOSEINC)
C
C  DSCALE is the distance scaling for decisions about whether to connect a given
C  pair of local minima. PSCALE is the scaling for the difference in GPFOLD.
C  IF (EXP(-(DISTANCE-DSCALE)/DSCALE).LT.RANDOM) is used in cycle.
C
      ELSE IF (WORD.EQ.'DSCALE') THEN
         CALL READF(DSCALE)
C
C  Setting DUMMYRUN to true means that no OPTIM jobs are actually submitted.
C
      ELSE IF (WORD.EQ.'DUMMYRUN') THEN
         DUMMYRUNT=.TRUE.
C
C  Setting DUMMYTST to true means that we create phoney entries in ts.data
C  corresponding to nearest-neighbour minima. For use with BHINTERP and BISECTT.
C
      ELSE IF (WORD.EQ.'DUMMYTS') THEN
         DUMMYTST=.TRUE.
C
C  DUMPGROUPS specifies that the groups of potential energy minima and transition
C  states obtained from subroutine regroupfree2 should be saved in files
C  minima_groups and ts_groups
C
      ELSE IF (WORD.EQ.'DUMPGROUPS') THEN
         DUMPGROUPST=.TRUE.
C
C  Energy difference  criterion for distinguishing stationary points
C  Can also be specified with ETOL
C
      ELSE IF (WORD.EQ.'EDIFFTOL') THEN
         CALL READF(EDIFFTOL)

!  Turns on the electric field   
      ELSE IF (WORD.EQ.'EFLD') THEN
         EFIELDT = .TRUE.
C
      ELSE IF (WORD.EQ.'ENERGY') THEN
         CALL READF(TOTALE)
         ENSEMBLE='E'
C
C  Energy difference  criterion for distinguishing stationary points
C
      ELSE IF (WORD.EQ.'ETOL') THEN
         CALL READF(EDIFFTOL)
C
C  Threshold for zero eigenvalues used to skip triples in GETALLPATHS
C
      ELSE IF (WORD.EQ.'EVCUT') THEN
         CALL READF(EVCUT)
C
C  OPTIM executable.
C
      ELSE IF (WORD.EQ.'EXEC') THEN
         CALL READA(EXEC)
C
C  GMIN executable.
C
      ELSE IF (WORD.EQ.'EXECGMIN') THEN
         CALL READA(EXECGMIN)
C
C  Write the coordinates of minimum WHICHMIN to file extractedmin and stop.
C
      ELSE IF (WORD.EQ.'EXTRACTMIN') THEN
         CALL READI(WHICHMIN)
         EXTRACTMINT=.TRUE.
C
C  Write the coordinates of minimum WHICHTS to file extractedts and stop.
C
      ELSE IF (WORD.EQ.'EXTRACTTS') THEN
         EXTRACTTST=.TRUE.
         CALL READI(WHICHTS)
C
C Choose connection pairs based on free energy barriers between regrouped minima
C and product minima.
C
      ELSE IF (WORD.EQ.'FREEPAIRS') THEN
         FREEPAIRT=.TRUE.
         CALL READF(REGROUPFREETHRESH)
         CALL READF(FREETHRESH)
         CALL READF(EINC)
         IF (EINC.LE.0.0D0) THEN
            PRINT '(2(A,G20.10))','keywords> WARNING EINC for FREEPAIRS reset from ',EINC,' to ',1.0D0
            EINC=1.0D0
         ENDIF
C
C  Frozen atoms (to adjust zeros when reading frequencies)
C  and needed for freezing atoms is tssearch.
C
      ELSE IF (WORD.EQ.'FREEZE') THEN
         IF (NATOMS.LE.0) THEN
            PRINT '(A,I6,A)','keywords> ERROR - NATOMS=',NATOMS,' NATOMS keyword must preceed FREEZE'
            STOP
         ENDIF
         IF (.NOT.ALLOCATED(FROZEN)) THEN
            ALLOCATE(FROZEN(NATOMS))
            DO J1=1,NATOMS
               FROZEN(J1)=.FALSE.
            ENDDO
         ENDIF
         FREEZE=.TRUE.
         DO J1=1,NITEMS-1
            NFREEZE=NFREEZE+1
            CALL READI(NDUMMY)
            FROZEN(NDUMMY)=.TRUE.
         ENDDO
         IF (PERMDIST) THEN
            NDUMMY=0
            DO J1=1,NPERMGROUP
               DO J2=1,NPERMSIZE(J1)
                  IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                     PRINT '(A,I8,A)',' keyword> ERROR atom ',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                     STOP
                  ENDIF
               ENDDO
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
         ENDIF       
C
C Friction coefficient for Berezhkovsii, Pollak, Zitserman formulation
C JCP, 97, 2422, 1992
C
      ELSE IF (WORD.EQ.'FRICTION') THEN
         FRICTIONT=.TRUE.
         IMFRQT=.TRUE.
         CALL READF(GAMMAFRICTION)
C
C  Distance criterion for distinguishing stationary points
C
      ELSE IF (WORD.EQ.'GEOMDIFFTOL') THEN
         CALL READF(GEOMDIFFTOL)
C
C  Specify graph transformation rate calculation. NCONNMIN is the connectivity
C  at which (and below) minima are removed. GTINT is the interval for performing
C  the analysis during a DPS calculation. IF GTINT=1 we do the analysis every cycle,
C  consisting of NCPU jobs, if GTINT=2 we do the analysis every other cycle. 
C  GTINT is zero by default, in which case we only do the analysis if NATTEMPTS=0.
C  GT is the original DJW implementation where we renormalise out all intervening
C  miinima. The resulting branching probablities and waiting times do not allow
C  for return to the starting minimum - they are like the rejection-free BKL approach.
C
      ELSE IF (WORD.EQ.'GT') THEN
         GTT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
         IF (NITEMS.GT.2) CALL READI(GTINT)
C
C  Semen Trygubenko, last modification Thu Mar 15 12:52:36 GMT 2007
C
C  A different implementation of GT method, as described in our GT JCP paper.
C  Input format: {SGT|DGT|SDGT} [DisconnectSources [AltPbb [Rescale [Normalise]]]]
C  Example input line: SDGT F T F F
C 
C  If set to true DisconnectSources instructs to disconnect sources, yielding rate
C  that theoretically corresponds to the one obtained stochastically with KMC. Disconnection
C  of the sources is currently implemented in DGT part of GT. For SGT runs with DisconnectSources=T
C  a call is made to DGT once all SGT work is done. DisconnectSources=T can be memory hungry.
C
C  If set to true, GT2Sparse instructs to perform GT analysis using
C  sparse-optimised algorithms and data structures.
C  
C  GT2Switch=T instructs the program to change the data structures and algorithms
C  from sparse to dense when the density requirement (can be adjusted by GT2RSWITCH) is
C  met.
C  
C  GT2AltPbb, when set to true, triggers the evaluation of Pbb sums using the algorithm
C  suggested by David to maintain precision. Alternative remedy is provided with
C  GT2Rescale which I have devised to stop the errors propagation. Both of these
C  should double the execution time, and are recognized by both
C  SGTDetachNode and DGTDetachNode routines. NB: Disconnect routine, another place where
C  roundoff errors can breed like rabbits, does not recognize either of these as of
C  Thu Mar 15 12:45:55 GMT 2007.
C  
C  GT2Normalise=T will instruct GT2input.f90 to normalize the branching probabilities when
C  obtained from PATHSAMPLE. This is purely for debugging purposes, as they always should be.
C 
      ELSE IF (word=='SGT'.or.word=='DGT'.or.word=='SDGT'.or.word=='GT2') then
         GT2T=.TRUE.
         if (word=='SGT') then
              GT2Sparse=.True.
              GT2Switch=.False.
         else if (word=='DGT') then
              GT2Sparse=.False.
              GT2Switch=.False.
         else if (word=='SDGT') then
              GT2Sparse=.True.
              GT2Switch=.True.
         else
              print '(A)','Please specify SGT, DGT or SDGT. GT2 keyword is now obsolete';stop
         endif
         PRINT *,'keywords> GT2Sparse=',GT2Sparse
         PRINT *,'keywords> GT2Switch=',GT2Switch
    
         IF (NITEMS>1) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.or.loyno=='F') then
                   READ(LOYNO,'(l1)') GT2DisconnectSources
                   PRINT '(1x,A,L5)','keywords> GT2DisconnectSources=',GT2DisconnectSources
              ELSE
                   PRINT '(1x,A,L5)','keywords> Invalid after "GT2" keyword'; stop
              ENDIF
         ENDIF
         IF (NITEMS>2) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.or.loyno=='F') then
                   READ(LOYNO,'(l1)') GT2AltPbb
                   PRINT '(1x,A,L5)','keywords> GT2AltPbb=',GT2AltPbb
              ELSE
                   PRINT '(1x,A,L5)','keywords> Invalid after "GT2" keyword'; stop
              ENDIF
         ENDIF
         IF (NITEMS>3) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.or.loyno=='F') then
                   READ(LOYNO,'(l1)') GT2Rescale
                   PRINT *,'keywords> GT2Rescale=',GT2Rescale
              ELSE
                   PRINT *,'keywords> Invalid after "GT2" keyword'; stop
              ENDIF
         ENDIF
         IF (NITEMS>4) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.or.loyno=='F') then
                   READ(LOYNO,'(l1)') GT2Normalise
                   PRINT *,'keywords> GT2Normalise=',GT2Normalise
              ELSE
                   PRINT *,'keywords> Invalid after "GT2" keyword'; stop
              ENDIF
         ENDIF
         IF (NITEMS>5) THEN
              print '(1x,a)','Input error: more than 4 entities following SGT, DGT or SDGT keyword!';stop
         ENDIF
C
C Switching ratio. Specifies when to switch from SGT to DGT.
C See GT2.f90 for details.
C
      ELSE IF (WORD=='GT2RSWITCH') then
         CALL READF(GT2RSWITCH)
C
C Is used to establish whether a node is a dead end or not.
C See GT2.f90 for details.
C
      ELSE IF (WORD=='GT2PTOL') then
         CALL READF(GT2PTOL)
C
C Prints the negative eigenvalue of each TS in ts.data as final (ninth) column
C
      ELSE IF (WORD=='IMFRQ') THEN
         IMFRQT=.TRUE.
C
C  Use constraint potential for interpolation as a connection metric (instead of distance).
C
      ELSE IF (WORD.EQ.'INTCONSTRAINT') THEN
         INTCONSTRAINTT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(INTCONSTRAINTTOL)
         IF (NITEMS.GT.2) CALL READF(INTCONSTRAINTDEL)
         IF (NITEMS.GT.3) CALL READF(INTCONSTRAINTREP)
         IF (NITEMS.GT.4) CALL READF(INTCONSTRAINREPCUT)
         IF (NITEMS.GT.5) CALL READI(INTCONSEP)
         IF (NITEMS.GT.6) CALL READI(INTREPSEP)
         INTERPCOSTFUNCTION=.TRUE.
C
C  Use interpolation potential for LJ.
C
      ELSE IF (WORD.EQ.'INTLJ') THEN
         INTLJT=.TRUE. 
         IF (NITEMS.GT.1) CALL READF(INTLJDEL)
         IF (NITEMS.GT.2) CALL READF(INTLJEPS)
         INTERPCOSTFUNCTION=.TRUE.
C
C  Inertia difference criterion - no longer used for distinguishing stationary points!
C
      ELSE IF (WORD.EQ.'ITOL') THEN
         CALL READF(IDIFFTOL)
C
C  Number of jobs to run per node for a distributed memory architecture.
C
      ELSE IF (WORD.EQ.'JOBSPERNODE') THEN
         CALL READI(JPN)
         CALL GETNODES(NCPU)
C
C  Semen Trygubenko, Thu Mar 15 16:21:41 GMT 2007
C  Minimum number of connections a minimum ought to have to be included in the
C  rate calculation or regrouping schemes. Can be set directly for some
C  keywords, but not others, so this keyword is provided just to make sure
C  it casn be set for all methods.
C
      ELSE IF (WORD.EQ.'NCONNMIN') THEN
         if (NITEMS.GT.1) then
              CALL READI(NCONNMIN)
         else
              print *, 'keywords> Usage: RateNConnMin <integer>'; stop
         endif
      ELSE IF (WORD.EQ.'KMC') THEN
         KMCT=.TRUE.
         NOPOINTS=.TRUE.
         IF (NITEMS.GT.1) CALL READF(NKMCCYCLES)
         IF (NITEMS.GT.2) CALL READF(PAIRTHRESH)
         IF (NITEMS.GT.3) CALL READI(NCONNMIN)
C
C  Specify new graph transformation rate calculation. NCONNMIN is the connectivity
C  at which (and below) minima are removed. 
C  NGT is the new implementation where we renormalise out all intervening
C  miinima but allow return to the starting state. This should give committor
C  probabilities if we remove all the I minima, unlike GT and GT2.
C
      ELSE IF (WORD.EQ.'NGT') THEN
         NGTT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
         IF (NITEMS.GT.2) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.OR.LOYNO=='F') then
                 READ(LOYNO,'(l1)') NGTDISCONNECTALL
                 IF (NGTDISCONNECTALL) PRINT '(A)','keywords> NGT will calculate kSS, kNSS and kKMC'
                 IF (.NOT.NGTDISCONNECTALL) PRINT '(A)','keywords> NGT will calculate kSS and kNSS but not kKMC'
              ELSE
                 PRINT '(1X,A,L5)','keywords> Invalid after "NGT" keyword'; STOP
              ENDIF
         ENDIF
         IF (NITEMS.GT.3) CALL READF(NGTSWITCH)
         IF (NITEMS.GT.4) CALL READI(NGTSIZE)
         PRINT '(A,F12.4)','keywords> NGT will switch to dense renormalisation scheme at threshold ',NGTSWITCH
         PRINT '(A,I6)','keywords> NGT maximum square matrix size in dense phase=',NGTSIZE
         IF (NITEMS.GT.5) CALL READF(NGTCRSWITCH) ! N.B. threshold in giga-bytes
         PRINT '(A,F12.4)','keywords> NGT will use compressed row storage scheme beyond threshold ',NGTCRSWITCH
C
C  Set KMC parameters: number of KMC runs for averages and value of PAIRTHRESH.
C  If the product of branching probabilities p12*p21 > PAIRTHRESH then we
C  renormalise this pair.
C
      ELSE IF (WORD.EQ.'KMCCOMMIT') THEN
         NOPOINTS=.TRUE.
         KMCCOMMITT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(NKMCCYCLES)
         IF (NITEMS.GT.2) CALL READF(MAXBREAK)
         IF (NITEMS.GT.3) CALL READF(PABCONV)
         IF (NITEMS.GT.4) CALL READF(OMEGA)
         IF (NITEMS.GT.5) CALL READF(PAIRTHRESH)
         IF (NITEMS.GT.6) CALL READI(NCONNMIN)
C
C  k-th shortest paths analysis for each B (or A) minimum.
C
      ELSE IF (WORD.EQ.'KSHORTESTPATHS') THEN
         KSHORTESTPATHST=.TRUE.
         CALL READI(NPATHS)
         CALL READI(NCONNMIN)
         IF (NITEMS.GT.3) THEN 
            CALL READA(LOYNO)
            IF (LOYNO == 'T') KSHORT_FULL_PRINTT = .TRUE.
         ENDIF
C
C  Whether to read extra curvatures from min.data.info files in DUMMYTS runs
C
      ELSE IF (WORD.EQ.'LOWESTFRQ') THEN
         LOWESTFRQT=.TRUE.

      ELSE IF (WORD.EQ.'MACHINE') THEN
         MACHINE=.TRUE.
C
C  MAXBARRIER requires both sides to be greater than MAXBARRIER to discard.
C
      ELSE IF (WORD.EQ.'MAXBARRIER') THEN
         CALL READF(MAXBARRIER)
C
C  The maximum number of constraints to use in the constrained potential.
C  The default is 3.
C
      ELSE IF (WORD.EQ.'MAXCON') THEN
         CALL READI(MAXCONUSE)
C
C  MAXDOWNBARRIER checks the downhill barrier in Dijkstra.
C
      ELSE IF (WORD.EQ.'MAXDOWNBARRIER') THEN
         CALL READF(MAXDOWNBARRIER)
C
C  TSTHRESH discards transition states above the specfied threshold. May be useful
C  for producing a better initial path and excluding CHARMM transition states with
C  silly energies from the database. MAXTSENERGY does the same as TSTHRESH keyword 
C  for compatability with OPTIM.
C
      ELSE IF (WORD.EQ.'MAXTSENERGY') THEN
         CALL READF(TSTHRESH)
C
C  Add the minima and transition= states in path.info.<PATHNAME> and
C  output.<PATHNAME> to an existing database. The end points are NOT assumed to
C  belong to the A and B sets.
C
      ELSE IF (WORD.EQ.'MERGEDB') THEN
         MERGEDBT=.TRUE.
         CALL READA(PATHNAME)

      ELSE IF (WORD.EQ.'MSSTOCK') THEN

         MSSTOCKT = .TRUE.
         CALL READI(NRBSITES)
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL DEFMULTSTOCK()

C
C  Number of atoms - essential unless RIGIDBODIES/RBAA keyword is present.
C
      ELSE IF (WORD.EQ.'NATOMS') THEN
         CALL READI(NATOMS)
      ELSE IF (WORD.EQ.'NINTS') THEN
         CALL READI(NINTS)
C
C  If NOFRQS is specified then frequencies are assumed to be absent from path.info files.
C
      ELSE IF (WORD.EQ.'NOFRQS') THEN
         NOFRQS=.TRUE.
C
C  If NOPOINTS is specified then setup should not try to read the min.points or ts.points
C  files. Should be the default for post-DPS database kinetics analysis such as KMC.
C
      ELSE IF (WORD.EQ.'NOPOINTS') THEN
         NOPOINTS=.TRUE.
C
C  Read in an order parameter threshold, ORDERPARAM, that tells us when we are in the other 
C  phase for a DOALLMIN run.
C
      ELSE IF (WORD.EQ.'ORDER') THEN
         CALL READF(ORDERPARAM)
C
C  Set the frequency with which the list of pairs for future searches is recreated.
C  Currently works for default search type based on Pfold difference.
C

      ELSE IF (WORD.EQ.'PAHA') THEN

         CALL READI(PAHID)

         IF (PAHID == 1) THEN
            NRBSITES = 12
         ENDIF

         PAHAT    = .TRUE.
         ALLOCATE(RBSITE(NRBSITES,3))

         IF (PAHID == 1) THEN
            CALL DEFBENZENE()
         ENDIF


      ELSE IF (WORD.EQ.'PAIRLIST') THEN
         CALL READI(NPAIRFRQ)

      ELSE IF (WORD.EQ.'PATCHYD') THEN
         PATCHYDT = .TRUE.
         CALL READI(NRBSITES)
         IF (NRBSITES .NE. 4) THEN
            PRINT *, 'NRBSITES has to be 4'
            STOP 
         ENDIF
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL DEFPATCHYD()
!
!  Whether to optimise the permutational isomers in assessing optimal alignment.
!  For PERMDIST all minimum distances will be minimsed with respect to the
!  specified permutations. For PERMISOMER we only check for permutational isomers
!  if the energy difference is below EDIFFTOL. This should save a great deal
!  of time for large systems containing many equivalent atoms, although the
!  distance won't be minimised with repect to permutations for inequivalent minima.
!
      ELSE IF ((WORD.EQ.'PERMDIST').OR.(WORD.EQ.'PERMISOMER')) THEN
         IF (NATOMS.LE.0) THEN
            PRINT '(A,I6,A)','keywords> ERROR - NATOMS=',NATOMS,' NATOMS keyword must preceed PERMDIST'
            STOP
         ENDIF
         IF (WORD.EQ.'PERMDIST') PERMDIST=.TRUE.
         IF (WORD.EQ.'PERMISOMER') PERMISOMER=.TRUE.
         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
         ALLOCATE(NPERMSIZE(3*NATOMS),PERMGROUP(3*NATOMS),NSETS(3*NATOMS),SETS(NATOMS,3))
!
!  The above dimensions were fixed at NATOMS because:
!  (a) Atoms were not allowed to appear in more than one group.
!  (b) The maximum number of pair exchanges associated with a group is three.
!
! However, for flexible water models we need to exchange all waters,
! and we can exchange H's attached to the same O. The dimension required
! becomes 3*NATOMS
!
         IF (PERMFILE) THEN
            OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
            READ(1,*) NPERMGROUP
            NDUMMY=1
            DO J1=1,NPERMGROUP
!
!  Sanity checks!
! 
               READ(1,*) NPERMSIZE(J1),NSETS(J1)
               IF (NSETS(J1).GT.3) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of secondary sets ',NSETS(J1),' is > 3'
                  STOP
               ENDIF
               IF (NDUMMY+NPERMSIZE(J1)-1.GT.3*NATOMS) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
!              READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SWAP1(PERMGROUP(J3),J2),J2=1,NSWAP(J1)),
!    &                                                            J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1)
               READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1),
     &                                                              J2=1,NSETS(J1))
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
            CLOSE(1)
         ELSE
            NPERMGROUP=1 ! ALL ATOMS CAN BE PERMUTED - DEFAULT
            NPERMSIZE(1)=NATOMS ! ALL ATOMS CAN BE PERMUTED - DEFAULT
            DO J1=1,NATOMS
               PERMGROUP(J1)=J1
            ENDDO
         ENDIF
         PRINT '(A,I6)','keywords> Number of groups of permutable atoms=',NPERMGROUP
         NDUMMY=1
         DO J1=1,NPERMGROUP
            PRINT '(A,3(I6,A))',' keyword> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &                                                 NSETS(J1),' additional atom sets:'
            WRITE(*,'(22I6)',ADVANCE='NO') PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1) 
            IF (NSETS(J1).GT.0) THEN
               WRITE(*,'(A)',ADVANCE='NO') ' with '
               DO J2=1,NSETS(J1)
                  DO J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1
                     WRITE(*,'(I6)',ADVANCE='NO') SETS(PERMGROUP(J3),J2)
                     IF (J3.LT.NDUMMY+NPERMSIZE(J1)-1) WRITE(*,'(A3)',ADVANCE='NO') ' / '
                  ENDDO 
                  IF (J2.LT.NSETS(J1)) WRITE(*,'(A3)',ADVANCE='NO') ' ; '
               ENDDO 
            ENDIF
            PRINT *,' '
            NDUMMY=NDUMMY+NPERMSIZE(J1)
         ENDDO
!
!  Another check.
!
         IF (NFREEZE.GT.0) THEN
            NDUMMY=0
            DO J1=1,NPERMGROUP
               DO J2=1,NPERMSIZE(J1)
                  IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                     PRINT '(A,I8,A)',' keyword> ERROR atom ',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                     STOP
                  ENDIF
               ENDDO
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
         ENDIF
C
C  Initial PERT parameter for geometry perturbations used in single-ended ts searches
C
      ELSE IF (WORD.EQ.'PERTURB') THEN
         CALL READF(PERTVALUE)
         PERTMIN=PERTVALUE/2.0D0
         PERTMAX=PERTVALUE*2.0D0
         IF (NITEMS.GT.2) CALL READF(PERTMIN)
         IF (NITEMS.GT.3) CALL READF(PERTMAX)
C
C  NPFOLD is the number of iterations per call to the global Pfold subroutine
C  and PFOLDINT is the frequency at which we call this GPFOLD calculation, 
C  i.e. PFOLDINT=1 means for every cycle, PFOLDINT=2 means every other cycle etc.
C
      ELSE IF (WORD.EQ.'PFOLD') THEN
         CALL READI(NPFOLD)
         CALL READI(PFOLDINT)
         CALL READF(OMEGA)
C
C  The value of the Planck constant in units of prevailing energy * seconds. 
C  This is needed for regrouped ts free energies. 
C  For example, for CHARMM we enter the temperature in kcal/mol (i.e. k_B * T),
C  so h needs to be in kcal/mol * seconds i.e. 9.546 * 10^(-14).
C
      ELSE IF (WORD.EQ.'PLANCK') THEN
         CALL READF(PLANCK)
      ELSE IF (WORD.EQ.'PSCALE') THEN
         CALL READF(PSCALE)
C
C  Try connecting minima that have the largest values of equilibrium occupation probability times
C  waiting time to minima on the fastest path.
C  Don't bother - UNTRAP with free energy sorting should do the job.
C
!     ELSE IF (WORD.EQ.'PTAU') THEN
!        PTAUT=.TRUE.
!        DIJKSTRAT=.TRUE.
      ELSE IF (WORD.EQ.'PULL') THEN
         PULLT=.TRUE.
         PRINT '(A)','keywords> Constant fulling force with 4 zero eigenvalues'
      ELSE IF (WORD.EQ.'RBAA') THEN
         CALL READI(NATOMS)
         NATOMS=NATOMS*2
         NTSITES  = (NATOMS/2)*NRBSITES
         IF (DBPTDT) NTSITES = (NATOMS/2-1)*NRBSITES + 4
         IF (NRBSITES == 0) THEN
            PRINT *, 'NRBSITES not yet defined'
            STOP
         ENDIF
         RBAAT = .TRUE.
      ELSE IF (WORD.EQ.'RBSYM') THEN
         RBSYMT=.TRUE.
         INQUIRE(FILE='rbsymops',EXIST=RBSYMTEST)
         IF (RBSYMTEST) THEN
            OPEN(UNIT=1,FILE='rbsymops',STATUS='OLD')
            READ(1,*) NRBGROUP
            ALLOCATE(RBOPS(4,NRBGROUP))
            READ(1,*) ((RBOPS(J1,J2),J1=1,4),J2=1,NRBGROUP)
            PRINT '(A,I6)',' keywords> number of symmetry operations for rigid body=',NRBGROUP
            DO J1=1,NRBGROUP
               PRINT '(A,I6)',' keywords> rigid-body symmetry operation', J1
               RBOPS(4,J1) = RBOPS(4,J1)*ATAN(1.D0)/45.D0
               PRINT '(3F20.10)',RBOPS(1:4,J1)
            ENDDO
         ELSE
            PRINT '(A)',' keywords> ERROR *** missing file rbsymops'
            STOP
         ENDIF

      ELSE IF (WORD.EQ.'READMIN') THEN
         READMINT=.TRUE.
         CALL READA(MINNAME)
C
C  Threshold for redefining the A/B/I sets on the basis of a PE superbasin analysis
C  at energy REGROUPTHRESH.
C
      ELSE IF (WORD.EQ.'REGROUP') THEN
         NOPOINTS=.TRUE.
         REGROUPT=.TRUE.
         CALL READF(REGROUPTHRESH)
C
C  Regroup on the basis of free energies (equivalent to intergroup rates).
C  Subsequent rate calculations use the new groups and free energies.
C
      ELSE IF (WORD.EQ.'REGROUPFREE') THEN
         NOPOINTS=.TRUE.
         REGROUPFREET=.TRUE.
         CALL READF(REGROUPFREETHRESH)
C
C  Regroup on the basis of free energies (equivalent to intergroup rates).
C  Subsequent rate calculations use pe stationary points, but the
C  A and B groups are expanded based on the free energy regrouping.
C
      ELSE IF (WORD.EQ.'REGROUPFREEAB') THEN
         NOPOINTS=.TRUE.
         REGROUPFREEABT=.TRUE.
         CALL READF(REGROUPFREETHRESH)
C
C  Threshold for regrouping and calculating free energies
C  on the basis of a superbasin analysis at PE REGROUPPETHRESH
C
      ELSE IF (WORD.EQ.'REGROUPPE') THEN
         NOPOINTS=.TRUE.
         REGROUPPET=.TRUE.
         CALL READF(REGROUPPETHRESH)
C
C  Threshold for regrouping and calculating free energies
C  on the basis of rates.
C
      ELSE IF (WORD.EQ.'REGROUPRATE') THEN
         NOPOINTS=.TRUE.
         REGROUPRATET=.TRUE.
         CALL READF(REGROUPRATETHRESH)
C
C  Remove stationary points specified in file min.remove and ts.remove
C
      ELSE IF (WORD.EQ.'REMOVESP') THEN
         REMOVESP=.TRUE.
C
C  Remove stationary points specified in file min.remove and ts.remove
C
      ELSE IF (WORD.EQ.'REMOVEUNCONNECTED') THEN
         REMOVEUNCONNECTEDT=.TRUE.
          IF (NITEMS.GT.1) CALL READA(UNCONNECTEDS)
C
C  Retain only stationary points specified in file min.retain and ts.retain
C
      ELSE IF (WORD.EQ.'RETAINSP') THEN
         RETAINSP=.TRUE.
C
C  Reweighting for reactant minima to allow stochastic sampling of reactant
C  in a GT calculation.
C
      ELSE IF (WORD.EQ.'REWEIGHT') THEN
         REWEIGHTT=.TRUE.
         CALL READI(NRWBINS)      ! number of bins in probability distribution
         CALL READI(NRWREACTANT)  ! number of reactant minima to choose
         CALL READA(RWENERGYFILE) ! name of file containing quench data
         IF (ALLOCATED(RWPROB)) DEALLOCATE(RWPROB)
         ALLOCATE(RWPROB(NRWBINS))
         OPEN(UNIT=1,FILE=TRIM(ADJUSTL(RWENERGYFILE)),STATUS='OLD')
         RWEMAX=-1.0D100
         RWEMIN=1.0D100
         NDUMMY=0
         DO
            READ(1,*,END=222) DUMMY
            IF (DUMMY.GT.RWEMAX) RWEMAX=DUMMY
            IF (DUMMY.LT.RWEMIN) RWEMIN=DUMMY
            NDUMMY=NDUMMY+1
         ENDDO
222      CONTINUE
         PRINT '(A,I8,2A)','keyword> ',NDUMMY,' energies read from file ',TRIM(ADJUSTL(RWENERGYFILE))
         RWBINWIDTH=(RWEMAX-RWEMIN)/NRWBINS
         REWIND(1)
         RWPROB(1:NRWBINS)=0.0D0
         DO J1=1,NDUMMY
            READ(1,*) DUMMY
            NDUMMY2=INT((DUMMY-RWEMIN-1.0D-10)/RWBINWIDTH) + 1
!           PRINT '(A,3G20.10,I6)','RWEMIN,RWEMAX,DUMMY,NDUMMY2=',RWEMIN,RWEMAX,DUMMY,NDUMMY2
            RWPROB(NDUMMY2)=RWPROB(NDUMMY2)+1
         ENDDO
         CLOSE(1)
         RWPROB(1:NRWBINS)=RWPROB(1:NRWBINS)/NDUMMY
         IF (DEBUG) THEN
            PRINT '(A)','keyword> bin energy ranges and probabilities:'
            DUMMY=0.0D0
            DO J1=1,NRWBINS
               PRINT '(I6,2F20.10,G20.10)',J1,RWEMIN+RWBINWIDTH*(J1-1),RWEMIN+RWBINWIDTH*J1,RWPROB(J1)
               DUMMY=DUMMY+RWPROB(J1)
            ENDDO
            PRINT '(A,G20.10)','keyword> sum of probabilities=',DUMMY
         ENDIF
C
C  Number of rigid bodies - essential unless NATOMS keyword is present.
C  All we should have to do is then set NATOMS equal to twice the number
C  of rigid bodies to get all the dimensions right.
C
      ELSE IF (WORD.EQ.'RIGIDBODIES') THEN
         CALL READI(NATOMS)
         NATOMS=NATOMS*2
C
C  Random number seed.
C
      ELSE IF (WORD.EQ.'SEED') THEN
         CALL READI(ISEED)
C
C  Try connecting closest minima that are further apart than a given
C  separation in terms of steps in the best path.
C
      ELSE IF (WORD.EQ.'SHORTCUT') THEN
         SHORTCUTT=.TRUE.
         DIJKSTRAT=.TRUE.
         CALL READI(MINSEP)
         IF (NITEMS.GT.2) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'BARRIER') BARRIERSHORT=.TRUE.
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'RATE') RATESHORT=.TRUE.

      ELSE IF (WORD.EQ.'SIS')THEN
!
! SIS epidemiological model
!
         SIST=.TRUE.
         CALL READI(SMAX)
         CALL READI(IMAX)
         CALL READI(POPSS)
         CALL READF(SISMU)
         CALL READF(SISKAPPA)
         CALL READF(SISBETA)
         ZSYM='C' ! just setting this randomly in case it shouldn't be undefined...
         NATOMS=1 ! Just setting this randomly in case it shouldn't be undefined...
         PRINT '(A,3I6,1X,3G15.10)','keywords> SIS parameters ',SMAX,IMAX,POPSS,SISMU,SISKAPPA,SISBETA

      ELSE IF (WORD.EQ.'ST') THEN

         STOCKAAT = .TRUE.
         NRBSITES = 1
         ALLOCATE(RBSITE(NRBSITES,3))
         RBSITE(1,:) = 0.D0
C
C  Node specification in nodes.info file in slurm format.
C
      ELSE IF (WORD.EQ.'SLURM') THEN
         SLURMT=.TRUE.
         CALL READI(JPN)
         CALL GETNODES(NCPU)
C
C  Make the initial min.A and min.B files using information in path.info.<PATHNAME> and
C  output.<PATHNAME>. This sets up a single A and a single B minimum, which are specified
C  by STARTMINA and STARTMINB. For path.info files in the DUMPALLPATHS format, these
C  will usually not be the first and last entries!
C
      ELSE IF (WORD.EQ.'STARTFROMPATH') THEN
         STARTFROMPATH=.TRUE.
         CALL READA(PATHNAME)
         CALL READI(STARTMINA)
         CALL READI(STARTMINB)
C
C  Similarly, STARTTRIPLES and ADDTRIPLES specify the triples format for
C  path.info files read using the STARTFROM and ADDPATH keywords.
C
      ELSE IF (WORD.EQ.'STARTTRIPLES') THEN
         STARTTRIPLES=.TRUE.
C
C  OPTIM system symbol, e.g. AX, LS, etc.
C
      ELSE IF (WORD.EQ.'SYSTEM') THEN
         CALL READA(ZSYM)
C
C  Number and mass of a tagged atom. 
C
      ELSE IF (WORD.EQ.'TAG') THEN
         IF (NATOMS.LT.1) THEN
            PRINT '(A)','keywords> ERROR - number of atoms must be set before TAG keyword in pathdata'
            STOP
         ENDIF
         IF (.NOT.ALLOCATED(TAGFAC)) THEN
            TAGT=.TRUE.
            ALLOCATE(TAGFAC(NATOMS),TAGNUM(NATOMS))
            TAGFAC(1:NATOMS)=1.0D0
            TAGNUM(1:NATOMS)=0
         ENDIF
         NTAG=NTAG+1
         CALL READI(TAGNUM(NTAG))
         CALL READF(TAGFAC(TAGNUM(NTAG)))
C
C  Canonical temperature.
C
      ELSE IF (WORD.EQ.'TEMPERATURE') THEN
         CALL READF(TEMPERATURE)
         ENSEMBLE='T'
C
C  NTFOLD is the maximum number of iterations for the Tfold subroutine.
C  Make it real to prevent overflow.
C
      ELSE IF (WORD.EQ.'TFOLD') THEN
         TFOLDT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
         IF (NITEMS.GT.2) CALL READF(NTFOLD)
         IF (NITEMS.GT.3) CALL READF(TFOLDTHRESH)
         IF (NITEMS.GT.4) CALL READF(TOMEGA)
C
C  If TRIPLES is .TRUE. we read path.info as min-sad-min triples, rather
C  than the traditional min-sad-min-sad-...
C  This keyword must be used if odata.connect uses DUMPALLPATHS rather than
C  DUMPPATH.
C
      ELSE IF (WORD.EQ.'TRIPLES') THEN
         TRIPLES=.TRUE.
C
C  TSTHRESH discards transition states above the specfied threshold. May be useful
C  for producing a better initial path and excluding CHARMM transition states with
C  silly energies from the database.
C
      ELSE IF (WORD.EQ.'TSTHRESH') THEN
         CALL READF(TSTHRESH)
C
C  Two-dimensional flatland.
C
      ELSE IF (WORD.EQ.'TWOD') THEN
         TWOD=.TRUE.
C
C jmc Set unres potential. NDIHE is the number of dihedral angles in the system, which
C is used in perturb and tssearch to perturb the system randomly
C
      ELSE IF (WORD.EQ.'UNRES') THEN
         UNRST=.TRUE.
         CALL READI(NDIHE)
         ZSYM='C'
C
C Choose connection pairs based on pe barriers to minima in the product set
C
      ELSE IF (WORD.EQ.'UNTRAP') THEN
         UNTRAPT=.TRUE.
         CALL READF(EINC)
         CALL READF(EUNTRAPTHRESH)
C
C Choose connection pairs based on the sequence read from file USEPAIRSFILE
C
      ELSE IF (WORD.EQ.'USEPAIRS') THEN
         USEPAIRST=.TRUE.
         CALL READA(USEPAIRSFILE)
         ! </kwd>
         ! }}}
!op226 End section {{{
!op226 <end>
      ELSE

         CALL REPORT('Unrecognized command '//WORD,.TRUE.)
      ENDIF

      CALL FLUSH(6,ISTAT)
      GOTO 190

      RETURN
      END
      !op226 </end>
      !op226 }}}
