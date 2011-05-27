! Fill in the Hessian matrix (two-sided numerical approach)
! {{{

        DO J = 1, 3*N
            QO(J) = QO(J) + DELTA

            include bln.cic.inc.f90
            include bln.grad.inc.f90
            FQ_PLUS=GRAD

            QO(J) = QO(J) - 2.0*DELTA

            include bln.cic.inc.f90
            include bln.grad.inc.f90
            FQ_MINUS=GRAD

            QO(J) = QO(J) + DELTA
    
            DO I = J, 3*N
                HESS(I,J) = (FQ_PLUS(I) -  FQ_MINUS(I))/(2.0*DELTA)
                HESS(J,I) = HESS(I,J)
            ENDDO
        ENDDO

! }}}

