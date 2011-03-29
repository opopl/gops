!***********************************************************************
!SUBROUTINES:                                                          *  
!            BUILD_RHO_SITE                                            *
!            BUILD_V_SITE	                                       *
!***********************************************************************
! FUNCTIONS:  RHO_POT,RHO_POT_D,RHO_POT_DD                             *
!             VPOT,VPOT_D,VPOT_DD				       *             
!             FEMBED,FEMBED_D,FEMBED_DD                                *
!             FCUT, FCUT_2,FCUT_DD                                     *
!             MFUNC,MFUNC_D, MFUNC_DD,                                 *
!             GFUNC,GFUNC_D,GFUNC_DD                                   *
!             HFUNC                                                    *
!             DELTA_DIRAC                                              *
!                                                                      *             
!****|******************************************************************|

!****|******************************************************************|
        DOUBLE PRECISION FUNCTION RHO_POT(IPOT,R)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION RC
        INCLUDE 'ACKLAND_SMA.H'  
        INCLUDE 'ACKLAND_MISHIN_CU.H' 
        COMMON /PARAM_CUT_OFF/RC	
	INTEGER IPOT		
	DOUBLE PRECISION  R
	
!DEBUG       	WRITE(*,*) '================================'
!DEBUG	WRITE(*,*) RC
!DEBIG	STOP 
		
        IF(IPOT.EQ.1) THEN   	
           RHO_POT=EXP(-2*Q*(R/RO-1.0D0))*FCUT(IPOT,R,RC,DELTA)	   
        ELSEIF(IPOT.EQ.2) THEN	
           RHO_POT=(RO/R)**(2*Q)*FCUT(IPOT,R,RC,DELTA)		   
        ELSEIF(IPOT.EQ.3) THEN	
            RHO_POT=EXP(-2*Q*(R/RO-1.0D0))*FCUT(IPOT,R,RC,DELTA)		   
	ELSEIF(IPOT.EQ.4) THEN			
	RHO_POT=GFUNC(R,A,BETA1,BETA2,R03,R04)*FCUT(IPOT,R,RC,H)	
	ELSEIF ((IPOT.EQ.5).OR.(IPOT.EQ.6)) THEN
	RHO_POT = FPSI(R)
	ELSEIF ((IPOT.EQ.12).OR.(IPOT.EQ.13)) THEN
	RHO_POT =  (   0.77718711248373D0 * (5.6D0-R)**4 &
	            -  0.48102928454986D0 * (5.6D0-R)**5 &  
		    +  0.14501312593993D0 * (5.6D0-R)**6 &
		    -  0.021292226813959D0* (5.6D0-R)**7 &
		    +  0.001220921762567D0* (5.6D0-R)**8) * HFUNC(5.6D0-R)
	ELSE
	WRITE(*,*) 'ERREUR DE IPOT'	
	
	ENDIF 
	   
        RETURN
        END    

!****|******************************************************************|	
!****|******************************************************************|	
       DOUBLE PRECISION FUNCTION RHO_POT_D(IPOT,R)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
        DOUBLE PRECISION RC     
        INCLUDE 'ACKLAND_SMA.H'  
        INCLUDE 'ACKLAND_MISHIN_CU.H' 
	COMMON /PARAM_CUT_OFF/RC
	DOUBLE PRECISION R		
	INTEGER IPOT	      


        IF(IPOT.EQ.1.) THEN     	   	     	        
      RHO_POT_D=-(2*Q)/RO*EXP(-2*Q*(R/RO-1.0D0))*FCUT(IPOT,R,RC,DELTA)+     &
               EXP(-2*Q*(R/RO-1.0D0))*FCUT_D(IPOT,R,RC,DELTA)
     
     
        ELSEIF(IPOT.EQ.2) THEN
    
        RHO_POT_D=-2*Q/R*(RO/R)**(2*Q)*FCUT(IPOT,R,RC,DELTA)+               &
               (RO/R)**(2*Q)*FCUT_D(IPOT,R,RC,DELTA)	
     
        ELSEIF(IPOT.EQ.3) THEN 
	   	     	        
      RHO_POT_D=-(2*Q)/RO*EXP(-2*Q*(R/RO-1.0D0))*FCUT(IPOT,R,RC,DELTA)+     &
               EXP(-2*Q*(R/RO-1.0D0))*FCUT_D(IPOT,R,RC,DELTA)	
	
        ELSEIF(IPOT.EQ.4) THEN 	
   			    
       RHO_POT_D=GFUNC_D(R,A,BETA1,BETA2,R03,R04)*FCUT(IPOT,R,RC,H)    &
     	        +GFUNC(R,A,BETA1,BETA2,R03,R04)*FCUT_D(IPOT,R,RC,H)   
	ELSEIF ((IPOT.EQ.5).OR.(IPOT.EQ.6)) THEN
	
	RHO_POT_D=FPSI_D(R)
	ELSEIF ((IPOT.EQ.12).OR.(IPOT.EQ.13)) THEN
	RHO_POT_D =  (- 0.77718711248373D0 * 4.D0* (5.6D0-R)**3 & 
	              + 0.48102928454986D0 * 5.D0* (5.6D0-R)**4 &
	              - 0.14501312593993D0 * 6.D0* (5.6D0-R)**5 &
		      + 0.021292226813959D0* 7.D0* (5.6D0-R)**6 &
		      - 0.001220921762567D0* 8.D0* (5.6D0-R)**7) * HFUNC(5.6D0-R)
	ELSE
	WRITE(*,*) 'ERREUR DE IPOT'
	ENDIF
		
		        
        RETURN
        END  
!****|******************************************************************|	
        DOUBLE PRECISION  FUNCTION RHO_POT_DD(IPOT,R)
	IMPLICIT  DOUBLE PRECISION(A-H,O-Z)  
	DOUBLE PRECISION RC      
        INCLUDE 'ACKLAND_SMA.H'  
        INCLUDE 'ACKLAND_MISHIN_CU.H' 
	COMMON /PARAM_CUT_OFF/RC	
	DOUBLE PRECISION R		
	INTEGER IPOT	      


       IF(IPOT.EQ.1) THEN   
       RHO_POT_DD=(2*Q/RO)**2*EXP(-2*Q*(R/RO-1.0D0))*FCUT(IPOT,R,RC,DELTA)        &
              -2.0D0*(2*Q/RO)*EXP(-2*Q*(R/RO-1.0D0))*FCUT_D(IPOT,R,RC,DELTA)      &
                             +EXP(-2*Q*(R/RO-1.0D0))*FCUT_DD(IPOT,R,RC,DELTA)      
	     
        ELSEIF(IPOT.EQ.2) THEN
        RHO_POT_DD=2*Q*(2*Q+1)/R**2*(RO/R)**(2*Q)*FCUT(IPOT,R,RC,DELTA)          &
                  -   2*(2*Q)/R*(RO/R)**(2*Q)*FCUT_D(IPOT,R,RC,DELTA)            &
                   + (RO/R)**(2*Q)*FCUT_DD(IPOT,R,RC,DELTA)
       
	
        ELSEIF(IPOT.EQ.3) THEN
       RHO_POT_DD=(2*Q/RO)**2*EXP(-2*Q*(R/RO-1.0D0))*FCUT(IPOT,R,RC,DELTA)         &
              -2.0D0*(2*Q/RO)*EXP(-2*Q*(R/RO-1.0D0))*FCUT_D(IPOT,R,RC,DELTA)       &
                             +EXP(-2*Q*(R/RO-1.0D0))*FCUT_DD(IPOT,R,RC,DELTA)           
	          	   	
        ELSEIF(IPOT.EQ.4) THEN 	
!****|******************************************************************|	   			    
       RHO_POT_DD=GFUNC_DD(R,A,BETA1,BETA2,R03,R04)*FCUT(IPOT,R,RC,H)     &    
             +2*GFUNC_D(R,A,BETA1,BETA2,R03,R04)*FCUT_D(IPOT,R,RC,H)      &
                + GFUNC(R,A,BETA1,BETA2,R03,R04)*FCUT_DD(IPOT,R,RC,H)  
!****|******************************************************************|	   
	ELSEIF ((IPOT.EQ.5).OR.(IPOT.EQ.6)) THEN
	RHO_POT_DD = FPSI_DD(R)
	
	ELSEIF ((IPOT.EQ.12).OR.(IPOT.EQ.13)) THEN
	RHO_POT_DD =  (  0.77718711248373D0 *12.D0* (5.6D0-R)**2 & 
	               - 0.48102928454986D0 *20.D0* (5.6D0-R)**3 &
	               + 0.14501312593993D0 *30.D0* (5.6D0-R)**4 &
		       - 0.021292226813959D0*42.D0* (5.6D0-R)**5 &
		       + 0.001220921762567D0*56.D0* (5.6D0-R)**6) * HFUNC(5.6D0-R)
	ELSE
	WRITE(*,*) 'ERREUR DE IPOT'
	ENDIF	        
        RETURN
        END  	
	
!****|******************************************************************|
        DOUBLE PRECISION FUNCTION VPOT(IPOT,R)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION RC
        INCLUDE 'ACKLAND_SMA.H'  
        INCLUDE 'ACKLAND_MISHIN_CU.H' 
	COMMON /PARAM_CUT_OFF/RC	
	INTEGER IPOT
	DOUBLE PRECISION MFUNC,R	
	
        IF(IPOT.EQ.1) THEN   	
           VPOT=EXP(-P*(R/RO-1.0D0))*FCUT(IPOT,R,RC,DELTA)
           VPOT=A_R*VPOT	   	    
        ELSEIF(IPOT.EQ.2) THEN	
           VPOT=(RO/R)**P*FCUT(IPOT,R,RC,DELTA)
           VPOT=A_R*VPOT	   		   
        ELSEIF(IPOT.EQ.3) THEN		
           VPOT=(RO/R)**P*FCUT(IPOT,R,RC,DELTA)	
           VPOT=A_R*VPOT	   			   	   
	ELSEIF(IPOT.EQ.4) THEN	
		
	VPOT=(E1*MFUNC(R,R01,ALPHA1)+E2*MFUNC(R,R02,ALPHA2)+DD)*      &
             FCUT(IPOT,R,RC,H)                                        &
            -S1*HFUNC(RS1-R)*(RS1-R)**4                               &
            -S2*HFUNC(RS2-R)*(RS2-R)**4                               &   
            -S3*HFUNC(RS3-R)*(RS3-R)**4	
	ELSEIF ((IPOT.EQ.5).OR.(IPOT.EQ.6)) THEN
	
	VPOT = 0.5D0*FVARPHI(R)
	
	ELSEIF (IPOT.EQ.12) THEN
	
	VPOT = EXP(  12.33339230761400D0	         &
	           - 10.84732196908600D0*R	         &
	           +  4.57335244245080D0*R**2            &
		   -  0.85266291445935D0*R**3)*            HFUNC(2.3D0-R) * HFUNC0(R-1.D0)  & 
		    
	         +(- 14.261501929757D0*   (3.5D0-R)**4   &
		   + 15.850036758176D0*   (3.5D0-R)**5   &
	           - 11.325102264291D0*   (3.5D0-R)**6   &
		   - 4.0971114831366D0*   (3.5D0-R)**7   &
	           + 3.6739378016909D0*   (3.5D0-R)**8  ) * HFUNC(3.5D0-R)*HFUNC(R-2.3D0)  &
		 +(  1.3066813393823D0*   (6.0D0-R)**4   &
		   - 0.60542710718094D0*  (6.0D0-R)**5   &
		   + 1.0055527194350D0 *  (6.0D0-R)**6   &
		   - 0.14918186777562D0*  (6.0D0-R)**7   &
		   + 0.032773112059590D0* (6.0D0-R)**8  ) * HFUNC(6.D0-R)*HFUNC(R-2.3D0)   &
		 +(  0.011433120304691D0* (7.6D0-R)**4   &
		   - 0.021982172508973D0* (7.6D0-R)**5   &
		   - 0.012542439692607D0* (7.6D0-R)**6   &
		   + 0.025062673874258D0* (7.6D0-R)**7   &
		   - 0.0075442887837418D0*(7.6D0-R)**8  ) * HFUNC(7.6D0-R)*HFUNC(R-2.3D0) 

	VPOT=0.5D0*VPOT
	ELSEIF (IPOT.EQ.13) THEN
        
	VPOT = EXP(  12.8822300381920D0    		     &
	           - 12.1838501578140D0*R  		     &
	           +  5.5998956281737D0*R**2		     &
		   -  1.0915156420318D0*R**3)               *  HFUNC(2.3D0-R)*HFUNC0(R-1.D0)  & 
		   
		 +(   8.4670497139946D0*     (3.5D0-R)**4    &
		   - 46.183472786003D0*      (3.5D0-R)**5    &
		   + 79.633499844770D0*      (3.5D0-R)**6    &
		   - 64.847634731465D0*      (3.5D0-R)**7    &
		   + 19.454623850774D0*      (3.5D0-R)**8  ) * HFUNC(3.5D0-R)*HFUNC(R-2.3D0)  & 
		   
		 +(-  0.097845860135187D0*   (6.0D0-R)**4    &
		   -  0.47537134413743D0*    (6.0D0-R)**5    &
		   -  0.00096806164225329D0* (6.0D0-R)**6    &
		   -  0.16355187497617D0*    (6.0D0-R)**7    &
		   -  0.00090914903435333D0* (6.0D0-R)**8  ) * HFUNC(6.D0-R)*HFUNC(R-2.3D0)   &   
		     
		 +(-  0.022038480751134D0*   (7.6D0-R)**4    &
		   -  0.060955465943384D0*   (7.6D0-R)**5    &
		   +  0.11573689045653D0*    (7.6D0-R)**6    &
		   -  0.062697675088029D0*   (7.6D0-R)**7    &
		   +  0.011273545085049D0*   (7.6D0-R)**8  ) * HFUNC(7.6D0-R)*HFUNC(R-2.3D0) 

        VPOT=0.5D0*VPOT
	ELSE
	WRITE(*,*) 'ERREUR DE IPOT'	
	
	ENDIF    
	   
        RETURN
        END    
!****|******************************************************************|
        DOUBLE PRECISION FUNCTION VPOT_D(IPOT,R)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
        INCLUDE 'ACKLAND_SMA.H'  
        INCLUDE 'ACKLAND_MISHIN_CU.H'	
	DOUBLE PRECISION RC
	COMMON /PARAM_CUT_OFF/RC	 
	DOUBLE PRECISION MFUNC,MFUNC_D,R		
	INTEGER IPOT	      


        IF(IPOT.EQ.1) THEN
      	   	     	        
        VPOT_D=-P/RO*DEXP(-P*(R/RO-1.0D0))*FCUT(IPOT,R,RC,DELTA)+    &
               DEXP(-P*(R/RO-1.0D0))*FCUT_D(IPOT,R,RC,DELTA)
        VPOT_D=A_R*VPOT_D     
     
        ELSEIF(IPOT.EQ.2) THEN
    
	VPOT_D=-P/R*(RO/R)**P*FCUT(IPOT,R,RC,DELTA)+                 &
               (RO/R)**P*FCUT_D(IPOT,R,RC,DELTA)
        VPOT_D=A_R*VPOT_D     
        ELSEIF(IPOT.EQ.3) THEN
    
	VPOT_D=-P/R*(RO/R)**P*FCUT(IPOT,R,RC,DELTA)+                 &
               (RO/R)**P*FCUT_D(IPOT,R,RC,DELTA)        	
        VPOT_D=A_R*VPOT_D    
        ELSEIF(IPOT.EQ.4) THEN 

      VPOT_D=(E1*MFUNC_D(R,R01,ALPHA1)+E2*MFUNC_D(R,R02,ALPHA2))*    &
            FCUT(IPOT,R,RC,H)                                        &
            +4*S1*HFUNC(RS1-R)*(RS1-R)**3                            &
            +4*S2*HFUNC(RS2-R)*(RS2-R)**3                            &    
            +4*S3*HFUNC(RS3-R)*(RS3-R)**3                            &
      +    (E1*MFUNC(R,R01,ALPHA1)+E2*MFUNC(R,R02,ALPHA2)+DD)*       &
            FCUT_D(IPOT,R,RC,H)     		 	   
	ELSEIF ((IPOT.EQ.5).OR.(IPOT.EQ.6)) THEN
	
	VPOT_D = 0.5D0*FVARPHI_D(R)
     
	ELSEIF (IPOT.EQ.12) THEN
	
	VPOT_D= EXP( 12.333392307614D0	   -  10.847321969086D0*R	    &
	            +4.5733524424508D0*R**2 -  0.85266291445935D0*R**3)     &
		*( -  10.847321969086D0	   +   4.5733524424508D0*2.D0*R     &
		   -  0.85266291445935D0*3.D0*R**2)                         &   
		    * HFUNC(2.3D0-R) * HFUNC(R-1.D0)                        &  
	         +(  14.261501929757D0*   4.D0*  (3.5D0-R)**3   &
		   - 15.850036758176D0*   5.D0*  (3.5D0-R)**4   &
	           + 11.325102264291D0*   6.D0*  (3.5D0-R)**5   &
		   + 4.0971114831366D0*   7.D0*  (3.5D0-R)**6   &
	           - 3.6739378016909D0*   8.D0*  (3.5D0-R)**7  ) * HFUNC(3.5D0-R)*HFUNC(R-2.3D0)  &
		 +(- 1.3066813393823D0*   4.D0*  (6.0D0-R)**3   &
		   + 0.60542710718094D0*  5.D0*  (6.0D0-R)**4   &
		   - 1.0055527194350D0*   6.D0*  (6.0D0-R)**5   &
		   + 0.14918186777562D0*  7.D0*  (6.0D0-R)**6   &
		   - 0.032773112059590*   8.D0*  (6.0D0-R)**7  ) * HFUNC(6.D0-R)*HFUNC(R-2.3D0)	&
		 +(- 0.011433120304691D0* 4.D0*  (7.6D0-R)**3   &
		   + 0.021982172508973D0* 5.D0*  (7.6D0-R)**4   &
		   + 0.012542439692607D0* 6.D0*  (7.6D0-R)**5   &
		   - 0.025062673874258D0* 7.D0*  (7.6D0-R)**6   &
		   + 0.0075442887837418D0*8.D0*  (7.6D0-R)**7  ) * HFUNC(7.6D0-R)*HFUNC(R-2.3D0) 

	VPOT_D=0.5D0*VPOT_D
	ELSEIF (IPOT.EQ.13) THEN
        
	VPOT_D = EXP( 12.882230038192D0      -   12.183850157814*R               &
	           +5.5998956281737D0*R**2 -    1.0915156420318D0*R**3)        &
		 *(-12.183850157814        +    5.5998956281737D0*2.D0*R       &
		   -1.0915156420318D0*3.D0*R**2)                               &
		   * HFUNC(2.3D0-R) * HFUNC0(R-1.D0)                            &
		 +(- 8.4670497139946D0*     4.D0*  (3.5D0-R)**3    &
		   + 46.183472786003D0*     5.D0*  (3.5D0-R)**4    &
		   - 79.633499844770D0*     6.D0*  (3.5D0-R)**5    &
		   + 64.847634731465D0*     7.D0*  (3.5D0-R)**6    &
		   - 19.454623850774D0*     8.D0*  (3.5D0-R)**7  ) * HFUNC(3.5D0-R)*HFUNC(R-2.3D0)  & 
		 +(+ 0.097845860135187D0*   4.D0*  (6.0D0-R)**3    &
		   + 0.47537134413743D0*    5.D0*  (6.0D0-R)**4    &
		   + 0.00096806164225329D0* 6.D0*  (6.0D0-R)**5    &
		   + 0.16355187497617D0*    7.D0*  (6.0D0-R)**6    &
		   + 0.00090914903435333D0* 8.D0*  (6.0D0-R)**7  ) * HFUNC(6.D0-R)*HFUNC(R-2.3D0)   &	  
		 +(+ 0.022038480751134D0*   4.D0*  (7.6D0-R)**3    &
		   + 0.060955465943384D0*   5.D0*  (7.6D0-R)**4    &
		   - 0.11573689045653D0*    6.D0*  (7.6D0-R)**5    &
		   + 0.062697675088029D0*   7.D0*  (7.6D0-R)**6    &
		   - 0.011273545085049D0*   8.D0*  (7.6D0-R)**7  ) * HFUNC(7.6D0-R)*HFUNC(R-2.3D0) 
	VPOT_D=0.5D0*VPOT_D

	ELSE
	WRITE(*,*) 'ERREUR DE IPOT'
	ENDIF	                         
    
        RETURN
        END    
	
!****|******************************************************************|
        DOUBLE PRECISION FUNCTION VPOT_DD(IPOT,R)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
        INCLUDE 'ACKLAND_SMA.H'  
        INCLUDE 'ACKLAND_MISHIN_CU.H' 
	DOUBLE PRECISION RC
	COMMON /PARAM_CUT_OFF/RC	
	DOUBLE PRECISION MFUNC,MFUNC_D,MFUNC_DD,R		
	INTEGER IPOT	      


        IF(IPOT.EQ.1) THEN   
       VPOT_DD=(P/RO)**2*EXP(-P*(R/RO-1.0D0))*FCUT(IPOT,R,RC,DELTA)      &
          -2.0D0*(P/RO)*EXP(-P*(R/RO-1.0D0))*FCUT_D(IPOT,R,RC,DELTA)     &
              +EXP(-P*(R/RO-1.0D0))*FCUT_DD(IPOT,R,RC,DELTA)
         VPOT_DD=A_R*VPOT_DD   
	     
        ELSEIF(IPOT.EQ.2) THEN
        VPOT_DD=P*(P+1)/R**2*(RO/R)**P*FCUT(IPOT,R,RC,DELTA)          &
                      -2*P/R*(RO/R)**P*FCUT_D(IPOT,R,RC,DELTA)        &
                            +(RO/R)**P*FCUT_DD(IPOT,R,RC,DELTA)       
        VPOT_DD=A_R*VPOT_DD  
	
        ELSEIF(IPOT.EQ.3) THEN
       VPOT_DD=P*(P+1)/R**2*(RO/R)**P*FCUT(IPOT,R,RC,DELTA)           &
                      -2*P/R*(RO/R)**P*FCUT_D(IPOT,R,RC,DELTA)        &
                            +(RO/R)**P*FCUT_DD(IPOT,R,RC,DELTA)       
         VPOT_DD=A_R*VPOT_DD  
	          	  
        ELSEIF(IPOT.EQ.4) THEN 

!****|******************************************************************|	   			    
      VPOT_DD=(E1*MFUNC_DD(R,R01,ALPHA1)+E2*MFUNC_DD(R,R02,ALPHA2))*FCUT(IPOT,R,RC,H)       &
       +    2*(E1*MFUNC_D(R,R01,ALPHA1)+E2*MFUNC_D(R,R02,ALPHA2))*FCUT_D(IPOT,R,RC,H)       &
       +     (E1*MFUNC(R,R01,ALPHA1)+E2*MFUNC(R,R02,ALPHA2)+DD)*FCUT_DD(IPOT,R,RC,H)        &   
       -      12*S1*HFUNC(RS1-R)*(RS1-R)**2            					    &
       -      12*S2*HFUNC(RS2-R)*(RS2-R)**2       					    &
       -      12*S3*HFUNC(RS3-R)*(RS3-R)**2      					          
     		 
!****|******************************************************************|	   
 	ELSEIF ((IPOT.EQ.5).OR.(IPOT.EQ.6)) THEN
	
	VPOT_DD = 0.5D0*FVARPHI_DD(R)
	
	ELSE
	WRITE(*,*) 'ERREUR DE IPOT'
	ENDIF	                         
    
        RETURN
        END    		
	
	  
!****|******************************************************************|	
	DOUBLE PRECISION  FUNCTION FCUT(IPOT,R,RC,WIDTH)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
	DOUBLE PRECISION R,RC,WIDTH,TEST,X
	
	IF(IPOT.EQ.1.OR.IPOT.EQ.2.OR.IPOT.EQ.3) THEN
           TEST=(R-RC)/WIDTH
	   IF(TEST.LT.-100.D0)  THEN
           FCUT=1.0D0	
	   ELSEIF(TEST.GT.100.D0) THEN
           FCUT=0.0D0	
	   ELSE	
	   X=TEST	
           FCUT=1.0D0/(1.0D0+EXP(X))
	   ENDIF
	ELSEIF(IPOT.EQ.4) THEN	   
	    TEST=(R-RC)
	   IF(TEST.GT.0.D0)  THEN
           FCUT=0.0D0
	   ELSEIF(TEST.LE.0.D0) THEN
	   X=TEST/WIDTH	   
           FCUT= X**4/(1+X**4)	   	   
	   ENDIF
	ELSEIF ((IPOT.EQ.5).OR.(IPOT.EQ.6).OR.(IPOT.EQ.12).OR.(IPOT.EQ.13)) THEN
	FCUT=1.D0   		   
	ENDIF
	   			
        RETURN
        END   	
	
!****|******************************************************************|	
	DOUBLE PRECISION  FUNCTION FCUT_D(IPOT,R,RC,WIDTH)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
	DOUBLE PRECISION R,RC,WIDTH,TEST,X
	
	IF(IPOT.EQ.1.OR.IPOT.EQ.2.OR.IPOT.EQ.3) THEN	
          TEST=(R-RC)/WIDTH
	   IF(TEST.LT.-100.D0)  THEN
             FCUT_D=0.0D0	
	   ELSEIF(TEST.GT.100.D0) THEN
             FCUT_D=0.0D0		
	    ELSE
	     X=TEST		
             FCUT_D=-EXP(X)/(1+EXP(X))**2
             FCUT_D=FCUT_D/WIDTH          
	   ENDIF		
	ELSEIF(IPOT.EQ.4) THEN	
	    TEST=(R-RC)
	   IF(TEST.GT.0.D0)  THEN
	   FCUT_D=0.0D0	
	   ELSEIF(TEST.LE.0.D0) THEN
	   X=TEST/WIDTH
           FCUT_D= 4*X**3/(1+X**4)**2  	 
           FCUT_D=FCUT_D/WIDTH  	     
	   ENDIF		
	ELSEIF ((IPOT.EQ.5).OR.(IPOT.EQ.6).OR.(IPOT.EQ.12).OR.(IPOT.EQ.13)) THEN
	FCUT_D=1.D0   		   
	ENDIF	
		
        RETURN
        END   	
	   		
!****|******************************************************************|	
	DOUBLE PRECISION FUNCTION FCUT_DD(IPOT,R,RC,WIDTH)
	IMPLICIT  DOUBLE PRECISION(A-H,O-Z)	
	DOUBLE PRECISION R,RC,WIDTH,TEST,X
	
	IF(IPOT.EQ.1.OR.IPOT.EQ.2.OR.IPOT.EQ.3) THEN	
          TEST=(R-RC)/WIDTH
	   IF(TEST.LT.-100.D0)  THEN
             FCUT_DD=0.0D0	
	   ELSEIF(TEST.GT.100.D0) THEN
             FCUT_DD=0.0D0		
	    ELSE
	     X=TEST		
             FCUT_DD= -EXP(X)/( (1.0D0+EXP(X) ) )**2 +2.0D0*EXP(2*X)/(1.0D0+EXP(X) )**3    
           FCUT_DD=FCUT_DD/WIDTH**2     
	   ENDIF		
	ELSEIF(IPOT.EQ.4) THEN	
	    TEST=(R-RC)
	   IF(TEST.GT.0.D0)  THEN
	   FCUT_DD=0.0D0	
	   ELSEIF(TEST.LE.0.D0) THEN
	   X=TEST/WIDTH	   
           FCUT_DD=12*X**2/(1+X**4)**2-32*X**6/(1+X**4)**3
           FCUT_DD=FCUT_DD/WIDTH**2 	      	   
	   ENDIF		
	ELSEIF ((IPOT.EQ.5).OR.(IPOT.EQ.6).OR.(IPOT.EQ.12).OR.(IPOT.EQ.13)) THEN
	FCUT_DD=1.D0   		   
	
	ENDIF	
		
        RETURN
        END  	
		   	 	 	
!****|******************************************************************|
        DOUBLE PRECISION FUNCTION FEMBED(IPOT,X)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
        DOUBLE PRECISION X
	INTEGER IPOT
	DOUBLE PRECISION RC
	COMMON /PARAM_CUT_OFF/RC	
        INCLUDE 'ACKLAND_SMA.H'  
        INCLUDE 'ACKLAND_MISHIN_CU.H'
	INCLUDE 'ACKLAND_MENDELEV_FE.H'	

		
	  IF(IPOT.EQ.1.OR.IPOT.EQ.2.OR.IPOT.EQ.3) THEN	
             FEMBED=X**(ALPHA)
             FEMBED=-B_A*FEMBED	  
	        
	  ELSEIF(IPOT.EQ.4) THEN	     
	  
	  IF(X.LE.1) THEN	  	  
	   FEMBED=F0+0.5*F2*(X-1)**2+Q1*(X-1)**3+Q2*(X-1)**4 +Q3*(X-1)**5 +Q4*(X-1)**6 
        
	  ELSEIF(X.GT.1) THEN	  	  	
	   FEMBED=(F0+0.5*F2*(X-1)**2+Q1*(X-1)**3+QQ1*(X-1)**4)/(1+QQ2*(X-1)**3)    	   
	  ENDIF	   	  
	  	  
	  ELSEIF(IPOT.EQ.5) THEN
	  
	  FEMBED = -DSQRT(X) + APHI*X*X + APHI2*X*X*X*X  
	  
	  ELSEIF(IPOT.EQ.6) THEN
	   IF (X.LT.1.D0) THEN
	     FEMBED = -APHI*DSQRT(X) - APHI2*(1-DSQRT(X))*DLOG(2.D0-X)/DLOG(2.D0) 
	     ELSE
	     FEMBED = -APHI*DSQRT(X) 
	   ENDIF

	  ELSEIF(IPOT.EQ.12) THEN
	   FEMBED = -     DSQRT(X)                                         &
	            - 1.9162462126235D0*1.D-7*(X-60.D0)**4*HFUNC(X-60.D0) &
	            + 4.6418727035037D0*1.D-7*(X-70.D0)**4*HFUNC(X-70.D0) &
	            + 6.6448294272955D0*1.D-7*(X-80.D0)**4*HFUNC(X-80.D0) &
		    - 2.0680252960229D0*1.D-6*(X-85.D0)**4*HFUNC(X-85.D0) &
		    + 1.1387131464983D0*1.D-6*(X-90.D0)**4*HFUNC(X-90.D0)
	  ELSEIF(IPOT.EQ.13) THEN
	   FEMBED = - DSQRT(X)                                         &
	            + 3.2283012597866D0*1.D-7*(X-60.D0)**4*HFUNC(X-60.D0) &		
		    - 1.1552813894483D0*1.D-6*(X-70.D0)**4*HFUNC(X-70.D0) &
		    + 2.3747280268355D0*1.D-6*(X-80.D0)**4*HFUNC(X-80.D0) &
		    - 2.0379550826523D0*1.D-6*(X-85.D0)**4*HFUNC(X-85.D0) & 
		    + 4.9758343293936D0*1.D-7*(X-90.D0)**4*HFUNC(X-90.D0)


	 END IF 
	RETURN
        END
!****|******************************************************************|
        DOUBLE PRECISION FUNCTION FEMBED_D(IPOT,X)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
	DOUBLE PRECISION X
	DOUBLE PRECISION RC
	COMMON /PARAM_CUT_OFF/RC	
        INCLUDE 'ACKLAND_SMA.H'  
        INCLUDE 'ACKLAND_MISHIN_CU.H'
	INCLUDE 'ACKLAND_MENDELEV_FE.H'
	
          
	  IF(IPOT.EQ.1.OR.IPOT.EQ.2.OR.IPOT.EQ.3) THEN	
             FEMBED_D=ALPHA*X**(ALPHA-1.0D0)
             FEMBED_D=-B_A*FEMBED_D	     	  
	  ELSEIF(IPOT.EQ.4) THEN	     	  
	  IF(X.LE.1) THEN	  	  
	   FEMBED_D=F2*(X-1)+3*Q1*(X-1)**2+4*Q2*(X-1)**3+5*Q3*(X-1)**4+6*Q4*(X-1)**5 
        
	  ELSEIF(X.GT.1) THEN	  	  	
	 FEMBED_D=(F2*(X-1)+3*Q1*(X-1)**2+4*QQ1*(X-1)**3)/(1+QQ2*(X-1)**3)       &
                  -3*QQ2*(X-1)**2*(F0+0.5*F2*(X-1)**2+Q1*(X-1)**3+QQ1*(X-1)**4)  &
                  /(1+QQ2*(X-1)**3)**2      
          
	      	   
	   ENDIF	   	  	  

	  ELSEIF (IPOT.EQ.5) THEN
	  
	  FEMBED_D = -0.5D0/DSQRT(X) + 2.D0*APHI*X  &
	              + 4.D0*APHI2*X*X*X


	  ELSEIF(IPOT.EQ.6) THEN
	  IF (X.LT.1.D0) THEN
	    FEMBED_D = -APHI*0.5D0/DSQRT(X)                   &
	               + APHI2*(0.5D0*DLOG(2.D0-X)/DSQRT(X)  &
	               + (1-DSQRT(X))/(2.D0-X))/DLOG(2.D0) 
	   ELSE
	    FEMBED_D = -0.5D0*APHI/DSQRT(X)     
          END IF 	    
	    
	  ELSEIF(IPOT.EQ.12) THEN
	   FEMBED_D = - 0.5D0 / DSQRT(X)                                 &
	            + 4.D0*(                                           &
		    - 1.9162462126235D0*1.D-7*(X-60.D0)**3*HFUNC(X-60.D0) &
	            + 4.6418727035037D0*1.D-7*(X-70.D0)**3*HFUNC(X-70.D0) &
	            + 6.6448294272955D0*1.D-7*(X-80.D0)**3*HFUNC(X-80.D0) &
		    - 2.0680252960229D0*1.D-6*(X-85.D0)**3*HFUNC(X-85.D0) &
		    + 1.1387131464983D0*1.D-6*(X-90.D0)**3*HFUNC(X-90.D0) &
		           )
	  ELSEIF(IPOT.EQ.13) THEN
	   FEMBED_D = - 0.5D0 / DSQRT(X)                                 &
	            + 4.D0*(                                           &
	            + 3.2283012597866D0*1.D-7*(X-60.D0)**3*HFUNC(X-60.D0) &		
		    - 1.1552813894483D0*1.D-6*(X-70.D0)**3*HFUNC(X-70.D0) &
		    + 2.3747280268355D0*1.D-6*(X-80.D0)**3*HFUNC(X-80.D0) &
		    - 2.0379550826523D0*1.D-6*(X-85.D0)**3*HFUNC(X-85.D0) & 
		    + 4.9758343293936D0*1.D-7*(X-90.D0)**3*HFUNC(X-90.D0) &
		           )
        
	
	ENDIF

        RETURN
        END
	
!****|******************************************************************|	
        DOUBLE PRECISION FUNCTION FEMBED_DD(IPOT,X)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)	
        INCLUDE 'ACKLAND_SMA.H'  
        INCLUDE 'ACKLAND_MISHIN_CU.H'
	DOUBLE PRECISION RC
	COMMON /PARAM_CUT_OFF/RC	
	DOUBLE PRECISION X,DENOM
	INCLUDE 'ACKLAND_MENDELEV_FE.H'
          
	  IF(IPOT.EQ.1.OR.IPOT.EQ.2.OR.IPOT.EQ.3) THEN	
             FEMBED_DD=ALPHA*(ALPHA-1)*X**(ALPHA-2.0D0)
             FEMBED_DD=-B_A*FEMBED_DD	     	  
	  ELSEIF(IPOT.EQ.4) THEN	     	  
	  IF(X<1) THEN
	   FEMBED_DD=F2+6*Q1*(X-1)+12*Q2*(X-1)**2+20*Q3*(X-1)**3+30*Q4*(X-1)**4 
        
	  ELSEIF(X>1) THEN
	 DENOM=1+QQ2*(X-1)**3
	 	  	  	  	
	 FEMBED_DD=(F2+6*Q1*(X-1)+12*QQ1*(X-1)**2)/DENOM                              &
         -6*QQ2*(X-1)**2*(F2*(X-1)+3*Q1*(X-1)**2+4*QQ1*(X-1)**3)/DENOM**2             &
         -6*QQ2*(X-1)*(F0+0.5*F2*(X-1)**2+Q1*(X-1)**3+QQ1*(X-1)**4)/DENOM**2          &
          +18*QQ2**2*(X-1)**4*(F0+0.5*F2*(X-1)**2+Q1*(X-1)**3+QQ1*(X-1)**4)/DENOM**3  
                        	   
	   ENDIF	   	  	  
        	    
	  ELSEIF (IPOT.EQ.5) THEN
	    FEMBED_DD = 0.25D0/DSQRT(X**3) + 2.D0*APHI  &
	              + 12.D0*APHI2*X*X
		      
          ELSEIF (IPOT.EQ.6) THEN
	   IF (X.LT.1.D0) THEN
	     FEMBED_DD = 0.25D0*APHI/DSQRT(X**3)  + APHI2/DLOG(2.D0)*   &
	             (                                              &
		       -0.25D0*DLOG(2.D0-X)/DSQRT(X**3)             &
		       -1.0/((2.D0-X)*DSQRT(X))                     &
		       +(1.0-DSQRT(X))/(2.D0-X)**2                  &
		      ) 
	     ELSE
	      FEMBED_DD = 0.25D0*APHI/DSQRT(X**3)     
            END IF 
	 
        ENDIF

        RETURN
        END
!****|******************************************************************|
        DOUBLE PRECISION FUNCTION DELTA_DIRAC(I,J)
        INTEGER I,J

         DELTA_DIRAC=0.0D0
        IF(I.EQ.J) THEN
         DELTA_DIRAC=1.0D0
        ENDIF
        RETURN
        END
!****|******************************************************************|	
	DOUBLE PRECISION FUNCTION HFUNC(X)
	DOUBLE PRECISION X
	
	IF(X.LT.0.0) THEN
	  HFUNC=0.0D0
	ELSEIF(X.GE.0.0) THEN	  
	  HFUNC=1.0D0	 
	 ENDIF
        RETURN
        END	 
!****|******************************************************************|	
	DOUBLE PRECISION FUNCTION HFUNC0(X)
	DOUBLE PRECISION X
	
	IF(X.LE.0.0) THEN
	  HFUNC0=0.0D0
	ELSEIF(X.GT.0.0) THEN	  
	  HFUNC0=1.0D0	 
	 ENDIF
        RETURN
        END	 
!****|******************************************************************|
	DOUBLE PRECISION FUNCTION MFUNC(R,R0,ALPHA)
	DOUBLE PRECISION R,R0,ALPHA
	
	MFUNC=EXP(-2*ALPHA*(R-R0))-2*EXP(-ALPHA*(R-R0))
		
        RETURN
        END	 
!****|******************************************************************|
	DOUBLE PRECISION FUNCTION MFUNC_D(R,R0,ALPHA)
	DOUBLE PRECISION R,R0,ALPHA
	
       MFUNC_D=-2*ALPHA*EXP(-2*ALPHA*(R-R0))+2*ALPHA*EXP(-ALPHA*(R-R0))
		
        RETURN
        END
!****|******************************************************************|	
        DOUBLE PRECISION FUNCTION MFUNC_DD(R,R0,ALPHA)
	DOUBLE PRECISION R,R0,ALPHA
	
         MFUNC_DD=+4*ALPHA**2*EXP(-2*ALPHA*(R-R0))-2*ALPHA**2*EXP(-ALPHA*(R-R0))
		
        RETURN
        END		
!****|******************************************************************|
	DOUBLE PRECISION FUNCTION GFUNC(X,A,B1,B2,X1,X2)
	DOUBLE PRECISION X,A,B1,B2,X1,X2
	
	GFUNC=A*EXP(-B1*(X-X1)**2)+EXP(-B2*(X-X2))
		
        RETURN
        END	 		
!****|******************************************************************|
	DOUBLE PRECISION FUNCTION GFUNC_D(X,A,B1,B2,X1,X2)
	DOUBLE PRECISION X,A,B1,B2,X1,X2
	
	GFUNC_D=-2*B1*(X-X1)*A*EXP(-B1*(X-X1)**2)-B2*EXP(-B2*(X-X2))
		
        RETURN
        END	
!****|******************************************************************|
	DOUBLE PRECISION FUNCTION GFUNC_DD(X,A,B1,B2,X1,X2)
	DOUBLE PRECISION X,A,B1,B2,X1,X2
	
	GFUNC_DD=-2*A*B1*EXP(-B1*(X-X1)**2)+4*A*B1**2*(X-X1)**2*EXP(-B1*(X-X1)**2)   &
                        +B2**2*EXP(-B2*(X-X2)) 		
        RETURN
        END	
	 		
!****|******************************************************************|
!****|******************************************************************|
!****|******************************************************************|


!****|******************************************************************|
!****|******************************************************************|
       DOUBLE PRECISION FUNCTION FPSI(X)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
       INCLUDE 'ACKLAND_MENDELEV_FE.H'
       
       TEMP = 0.D0
       DO I=1,NPSI
         TEMP = TEMP + AP(I)*HFUNC(RP(I)-X)*(RP(I)-X)**3
       END DO
       FPSI = TEMP
       RETURN
       END  
!****|******************************************************************|
       DOUBLE PRECISION FUNCTION FPSI_D(X)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
       INCLUDE 'ACKLAND_MENDELEV_FE.H'
       
       TEMP = 0.D0
       DO I=1,NPSI
         TEMP = TEMP - 3.D0*AP(I)*HFUNC(RP(I)-X)*(RP(I)-X)**2
       END DO
       FPSI_D = TEMP
       RETURN
       END  
!****|******************************************************************|
       DOUBLE PRECISION FUNCTION FPSI_DD(X)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
       INCLUDE 'ACKLAND_MENDELEV_FE.H'
       
       TEMP = 0.D0
       DO I=1,NPSI
         TEMP = TEMP + 6.D0*AP(I)*HFUNC(RP(I)-X)*(RP(I)-X)
       END DO
       FPSI_DD = TEMP
       RETURN
       END  
!****|******************************************************************|
!****|******************************************************************|
       DOUBLE PRECISION FUNCTION FVARPHI (X)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'ACKLAND_MENDELEV_FE.H'
        
	 THREEI   =1.D0/3.D0      
        ZNFE2   = ZNFE*ZNFE
       AU_TO_EV = HART
        AU_TO_A = ABOHR
         TEMP   = 0.D0

       RS = 0.88534D0*ABOHR*ZNFE**(-THREEI)/DSQRT(2.D0)
       RX = X/RS
       
        IF (X.LT.R1) THEN 
          FVARPHI = ZNFE2 * FPHI(RX) * AU_TO_A* AU_TO_EV/X 
        ELSE IF ((X.GE.R1).AND.(X.LT.R2)) THEN 
          FVARPHI = DEXP (  BFE0 + BFE1*X + BFE2*X*X + BFE3*X*X*X )
        ELSE IF (X.GE.R2) THEN
          DO I=1,NVARPHI
             TEMP = TEMP + AF(I)*HFUNC(RF(I)-X)*(RF(I)-X)**3
          END DO
          FVARPHI = TEMP
        END IF              
        
       RETURN
       END 
!****|******************************************************************|
       DOUBLE PRECISION FUNCTION FVARPHI_D (X)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'ACKLAND_MENDELEV_FE.H'
         THREEI   = 1.D0/3.D0
        ZNFE2   = ZNFE*ZNFE
       AU_TO_EV = HART
        AU_TO_A = ABOHR
         TEMP   = 0.D0

       RS = 0.88534D0*ABOHR*ZNFE**(-THREEI)/DSQRT(2.D0)
       RX = X/RS

        IF (X.LT.R1) THEN 
          FVARPHI_D =  - ZNFE2*FPHI(RX)*AU_TO_EV*AU_TO_A/(X*X) &
                 + ZNFE2 * FPHI_D(RX) * AU_TO_EV*AU_TO_A / (X*RS)     
        ELSE IF ((X.GE.R1).AND.(X.LT.R2)) THEN 
          FVARPHI_D = (BFE1 + 2.0D0*BFE2*X + 3.0D0*BFE3*X*X) &
                     *DEXP (BFE0 + BFE1*X + BFE2*X*X + BFE3*X*X*X)
        ELSE IF (X.GE.R2) THEN
          DO I=1,NVARPHI
             TEMP = TEMP - 3.D0*AF(I)*HFUNC(RF(I)-X)*(RF(I)-X)**2
          END DO
          FVARPHI_D = TEMP
        END IF              
        
       RETURN
       END
!****|******************************************************************|
       DOUBLE PRECISION FUNCTION FVARPHI_DD (X)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'ACKLAND_MENDELEV_FE.H'
         THREEI   = 1.D0/3.D0
        ZNFE2   = ZNFE*ZNFE
       AU_TO_EV = HART
        AU_TO_A = ABOHR
         TEMP   = 0.D0

       RS = 0.88534D0*ABOHR*ZNFE**(-THREEI)/DSQRT(2.D0)
       RX = X/RS

        IF (X.LT.R1) THEN 
	FVARPHI_DD = 0.5D0*ZNFE2*FPHI(RX)*AU_TO_EV*AU_TO_A/(X*X*X) &
	             - ZNFE2*FPHI_D(RX)*AU_TO_EV*AU_TO_A/(X*X*RS)  &
		     - ZNFE2 * FPHI_D(RX) * AU_TO_EV*AU_TO_A / (X*X*RS) &
		     + ZNFE2 * FPHI_DD(RX) * AU_TO_EV*AU_TO_A / (X*RS*RS)	 		 
        ELSE IF ((X.GE.R1).AND.(X.LT.R2)) THEN 
          FVARPHI_DD = (2.0D0*BFE2 + 6.0D0*BFE3*X) &
                     *DEXP (BFE0 + BFE1*X + BFE2*X*X + BFE3*X*X*X)&
		     + &
                    (BFE1 + 2.0D0*BFE2*X + 3.0D0*BFE3*X*X)**2 &
                     *DEXP (BFE0 + BFE1*X + BFE2*X*X + BFE3*X*X*X)		     
        ELSE IF (X.GE.R2) THEN
          DO I=1,NVARPHI
             TEMP = TEMP + 6.D0*AF(I)*HFUNC(RF(I)-X)*(RF(I)-X)
          END DO
          FVARPHI_DD = TEMP
        END IF              
        
       RETURN
       END
!****|******************************************************************|
!****|******************************************************************|
!****|******************************************************************|
       DOUBLE PRECISION FUNCTION FPHI(X)
       DOUBLE PRECISION X
  
        FPHI = 0.1818D0*DEXP(-3.2D0*X) & 
            +  0.5099D0*DEXP(-0.9423D0*X) &
            +  0.2802D0*DEXP(-0.4029D0*X) &
            +  0.02817*DEXP(-0.2016*X)

       RETURN
       END
!****|******************************************************************|
       DOUBLE PRECISION FUNCTION FPHI_D(X)
       DOUBLE PRECISION X
  
       FPHI_D = -  0.58176D0*DEXP(-3.2D0*X) &
                -  0.480479D0*DEXP(-0.9423D0*X) & 
                -  0.112893D0*DEXP(-0.4029D0*X) &
                -  0.00567907D0*DEXP(-0.2016*X)

       RETURN
       END

!****|******************************************************************|
       DOUBLE PRECISION FUNCTION FPHI_DD(X)
       DOUBLE PRECISION X
  
       FPHI_DD =    1.86163D0*DEXP(-3.2D0*X) &
                +  0.452755D0*DEXP(-0.9423D0*X) & 
                +  0.0454846D0*DEXP(-0.4029D0*X) &
                +  0.0011449D0*DEXP(-0.2016*X)

       RETURN
       END

!****|******************************************************************|
!****|******************************************************************|
!****|******************************************************************|
!****|******************************************************************|

      SUBROUTINE BUILD_RHO_SITE(IPOT,RHO_SITE,RN,NDIR,NAT_UP,NDIR_MAX)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DOUBLE PRECISION RN(NAT_UP,NDIR_MAX,3), RHO_SITE(NAT_UP),    &
                          RHO_TEMP,NORMR,R,RTEMP(3)
        INTEGER NDIR(NAT_UP),NAT_UP,NDIR_MAX,IPOT    

  
       DO I=1,NAT_UP  
          RHO_TEMP=0.0D0
           DO NI=1,NDIR(I)
	    RTEMP(:)=RN(I,NI,:)	   
            R=SQRT(DOT_PRODUCT(RTEMP,RTEMP)) 
            RHO_TEMP=RHO_TEMP+RHO_POT(IPOT,R)	    
           END DO

            RHO_SITE(I)=RHO_TEMP	
        END DO
 
       RETURN
       END
!****|******************************************************************|       
       SUBROUTINE BUILD_V_SITE(IPOT,V_SITE,RN,NDIR,NAT_UP,NDIR_MAX)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DOUBLE PRECISION RN(NAT_UP,NDIR_MAX,3),V_SITE(NAT_UP),   &
                         V_TEMP,NORMR,R,RTEMP(3)
        INTEGER NDIR(NAT_UP),NAT_UP,NAT_MAX,NDIR_MAX,IPOT

  
       DO I=1,NAT_UP  
          V_TEMP=0.0D0
           DO NI=1,NDIR(I)	   
	    RTEMP(:)=RN(I,NI,:)	   
            R=SQRT(DOT_PRODUCT(RTEMP,RTEMP))           	   
            V_TEMP=V_TEMP+VPOT(IPOT,R)	    
           END DO
              V_SITE(I)=V_TEMP	
        END DO
 
       RETURN
       END      

      

