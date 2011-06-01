
! write to LFH: Quench number, Energy, Steps, RMS, Markov E, t (elapsed time) {{{

WRITE(LFH,111)  ' Quench number,     NQ=',              NQ,                        	&
                ' Computed energy,   E=',               E,	&
                ' Iterations performed, ITERATIONS=',   ITERATIONS,	&
                ' RMS=',                      RMS,	&
                ' Markov E=',                 E,	&
                ' Time per quench, t=',                        TIME-TSTART
! }}}

