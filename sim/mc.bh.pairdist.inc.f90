
!     csw34> PAIRDIST prints the distance between specified atom pairs
!     after each quench. 

IF (PAIRDISTT) THEN
!     Write end of previous line as using ADVANCE="NO"
                WRITE(MYPUNIT,*) " "
!     Write current quench number
                WRITE(MYPUNIT,'(I10)',ADVANCE="NO") NQ
!     For each pair, assign ATOM1 and ATOM2 arrays containing coordinates
             DO PAIRCOUNTER=1,NPAIRS
                ATOM1(:)=COORDS(3*PAIRDIST(PAIRCOUNTER,1)-2:3*PAIRDIST(PAIRCOUNTER,1))
                ATOM2(:)=COORDS(3*PAIRDIST(PAIRCOUNTER,2)-2:3*PAIRDIST(PAIRCOUNTER,2))
!     Call PAIRDISTANCE with (x,y,z) for each atom
                WRITE(MYPUNIT,'(F10.4)',ADVANCE="NO") PAIRDISTANCE(ATOM1,ATOM2) 
             ENDDO
             FLUSH(MYPUNIT)
ENDIF
          

