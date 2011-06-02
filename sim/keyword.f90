
      SUBROUTINE KEYWORD(PROG)
      ! declarations {{{
      INTEGER PROG

      ! INPUT RELATED VARIABLES

      CHARACTER(LEN=100) :: BUFFER, LABEL
      INTEGER :: POS, LINE
      INTEGER, PARAMETER :: DATA_FH = 15
      INTEGER :: IOS = 0

      COMMON /BUFINF/ BUFFER,POS,LABEL,IOS
      ! }}}

      SELECTCASE(PROG)
        CASE(1) ! GMIN {{{

      CALL OPENF(DATA_FH,"<","data")

! ios<0 end of file;
! ios=0 
! ios>0 error 

      DO WHILE (IOS == 0)
        READ(DATA_FH, '(A)', IOSTAT=IOS) BUFFER
        IF (IOS == 0) THEN

           POS = SCAN(BUFFER, '    ')
           LABEL = BUFFER(1:POS)
           CALL PARSE(BUFFER,' ',ARGS,NARGS)
           BUFFER = BUFFER(POS+1:)
           LABEL=ARGS(1)
   
          SELECTCASE(LABEL)
              ! {{{
    	      CASE('  ','NOTE','COMMENT','\\','#')
    	      CASE('STOP')
    	         RETURN
                 ! A-L {{{
    	      CASE('ACCRAT')
    	         CALL READ(BUFFER, *, IOSTAT=IOS) ACCRAT
    	      CASE('FQMAX','TIGHTCONV')
    	         CALL READ(BUFFER, *, IOSTAT=IOS) FQMAX
              CASE('EDIFF')
    	         READ(BUFFER, *, IOSTAT=IOS) EDIFF
              CASE('G46')
                 G46=.TRUE.
                 BLNT=.TRUE.
                 ! }}}
                 ! M-S {{{
              CASE('MAXIT')
    	         READ(BUFFER, *, IOSTAT=IOS) MAXIT
              CASE('M_LBFGS')
    	         READ(BUFFER, *, IOSTAT=IOS) M_LBFGS
              CASE('NACCEPT')
    	         READ(BUFFER, *, IOSTAT=IOS) NACCEPT
              CASE('NSAVE')
    	         READ(BUFFER, *, IOSTAT=IOS) NSAVE
              CASE('P46')
                 P46=.TRUE.
                 BLNT=.TRUE.
              CASE('PULL')
                 PULLT=.TRUE.
    	         READ(BUFFER, *, IOSTAT=IOS) PATOM1,PATOM2,PFORCE
              CASE('RADIUS')
    	         READ(BUFFER, *, IOSTAT=IOS) RADIUS
    	      CASE('SQMAX','BASIN','SLOPPYCONV')
    	         CALL READ(BUFFER, *, IOSTAT=IOS) SQMAX
              CASE('STEPS')
    	         READ(BUFFER, *, IOSTAT=IOS) MCSTEPS,TFAC
                 ! }}}
              CASE('TARGET')
                TARGET=.TRUE.
                NTARGETS=NARGS-1
                ALLOCATE(TARGETS(NTARGETS))
    	        CALL READ(BUFFER, *, IOSTAT=IOS) TARGETS(1:NTARGETS)
              CASE('TRACKDATA')
                TRACKDATAT=.TRUE.
              CASE(DEFAULT)
    	         CALL REPORT('Unrecognized command '//BUFFER,.TRUE.)
    	         STOP
                 ! }}}
          ENDSELECT
       END IF
      END DO
! }}}
      ENDSELECT

      RETURN
      END

