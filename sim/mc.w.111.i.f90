
! write to LFH: Quench number, Energy, Steps, RMS, Markov E, t (elapsed time) {{{

WRITE(LFH,111)  ' Quench number,     NQ=',              NQ,                        	&
                ' Computed energy,   E=',               QE,	&
                ' ITERATIONS=',   ITERATIONS,	&
                ' RMS=',                      RMS,	&
                ' Markov E=',                 QE,	&
                ' Time per quench, t=',                        TIME-TSTART
! }}}

! write to LFH: Quench number, Energy, Steps, RMS, Markov E, t (elapsed time) {{{

WRITE(*,111)  ' Quench number,     NQ=',              NQ,                        	&
                ' Computed energy,   E=',               QE,	&
                ' ITERATIONS=',   ITERATIONS,	&
                ' RMS=',                      RMS,	&
                ' Markov E=',                 QE,	&
                ' Time per quench, t=',                        TIME-TSTART
! }}}

