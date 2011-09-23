
! Parameters for the non-bonded LJ interaction {{{

      AB=0.0D0
        DO I = 1, N
           DO J = 1, N
           IF ( ABS(I-J).GT.1) THEN

           IF (NTYPE(I) .EQ. 3 .OR. NTYPE(J) .EQ. 3) THEN
             AB(I,J,1) = 1.0*EPSILON 
             AB(I,J,2) = 0.0 
           ELSEIF (NTYPE(I) .EQ. 1 .AND. NTYPE(J) .EQ. 1)THEN
             AB(I,J,1) =  EPSILON
             IF (CONNECT(I,J)) THEN
                AB(I,J,2) = -EPSILON 
             ELSE
                AB(I,J,2) = 0.0D0
             ENDIF
           ELSE
             AB(I,J,1:2) = EPSILON*2.0/3.0 
           ENDIF
           ENDIF
   
           ENDDO
        ENDDO
! }}}
