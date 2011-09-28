
        IF (DEB) THEN 
		        NCALLMAX=100
		        IF (NCALL .LE. 0) THEN
		            NCALL=1
		          ELSE
		            NCALL=1+NCALL
		        ENDIF
		        IF (NCALL .GE. NCALLMAX) STOP
        ENDIF

