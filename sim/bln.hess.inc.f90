! Fill in the Hessian matrix
! {{{

        DO J = 1, 3*N
            QO(J) = QO(J) + DELTA

            CALL CALC_INT_COORDS(N,QO,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)
            CALL CALC_GRADIENT(N,QO,FQ_PLUS,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)

            QO(J) = QO(J) - 2.0*DELTA

            CALL CALC_INT_COORDS(N,QO,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)
            CALL CALC_GRADIENT(N,QO,FQ_MINUS,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)

            QO(J) = QO(J) + DELTA
    
            DO I = J, 3*N
                HESS(I,J) = (FQ_PLUS(I) -  FQ_MINUS(I))/(2.0*DELTA)
                HESS(J,I) = HESS(I,J)
            ENDDO
        ENDDO

! }}}

