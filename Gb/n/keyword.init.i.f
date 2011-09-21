
		! files {{{
		LFH=FH+1
		ENERGY_FH=FH+2
		MARKOV_FH=FH+3
		BEST_FH=FH+4
		PAIRDIST_FH=FH+5
		COORDS_FH=FH+6
		DATA_FH=FH+7
		RMSD_FH=FH+8
		! }}}

      NPCOUNT=0
      NPCALL=0
      NSEED=0
      NS=0
      NSSTOP=0
      NSAVE=5
      TFAC(:)=1.0D0
      ALLOCATE(FIXSTEP(1),FIXTEMP(1),FIXBOTH(1),TEMP(1),ACCRAT(1),STEP(1),ASTEP(1),OSTEP(1),NQ(1),EPREV(1),
     @         COORDS(3*NATOMS,1),COORDSO(3*NATOMS,1),VAT(NATOMS,1),VATO(NATOMS,1),
     @         NCORE(1))
      DO JP=1,1
         FIXSTEP(JP)=.FALSE.
         FIXTEMP(JP)=.FALSE.
         FIXBOTH(JP)=.FALSE.
         TEMP(JP)=0.3D0
         ACCRAT(JP)=0.5D0
         STEP(JP)=0.3D0
         ASTEP(JP)=0.3D0
         OSTEP(JP)=0.3D0
         NCORE(JP)=0
      ENDDO
      EDIFF=0.02D0
      DUMPT=.FALSE.
      TARGET=.FALSE.
      NTARGETS=0
      P46=.FALSE.
      G46=.FALSE.
      BLNT=.FALSE.
      DEBUG=.FALSE.
      SEEDT=.FALSE.


      CENT=.FALSE.
      CENTXY=.FALSE.
      CENTX=0.0D0
      CENTY=0.0D0
      CENTZ=0.0D0
      FIXCOM=.FALSE.

      PERIODIC=.FALSE.
      NRUNS=0
      RADIUS=0.0D0
      MAXIT=500
      MAXIT2=500
      FQMAX=1.0D-10
      SQMAX=1.0D-3
      NACCEPT=50
      NORESET=.FALSE.
      TRACKDATAT=.FALSE.
      MUPDATE=4
      DGUESS=0.1D0
      BFGS=.FALSE.
      LBFGST=.TRUE.
      RESTART=.FALSE.
      NRELAX=0
      NMSBSAVE=0
      NHSRESTART=0
      MAXBFGS=0.4D0

      EPSSPHERE=0.0D0
      MAXERISE=1.0D-10
      MAXEFALL=-HUGE(1.0D0)
      ARMA=0.4D0
      ARMB=0.4D0

      DUMPINT=1000 ! default is to dump a restart file every 1000 cycles of mc.f
      DUMPFILE=''
      RESTORET=.FALSE.
      INTEDUMPFILE=''

      PULLT=.FALSE.

      CHECKMARKOVT=.FALSE.


