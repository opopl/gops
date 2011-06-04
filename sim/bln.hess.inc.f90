! Fill in the Hessian matrix (two-sided numerical approach)
! {{{

        DO JH = 1, N
            DO KH=1,3
            R(JH,K) = R(JH,K) + DELTA

            include "bln.cic.inc.f90"
            include "bln.grad.inc.f90"
            GRAD_PLUS=G

            R(JH,K) = R(JH,KH) - 2.0*DELTA

            include "bln.cic.inc.f90"
            include "bln.grad.inc.f90"
            GRAD_MIN=G

            R(JH,KH) = R(JH,KH) + DELTA
    
            DO IH = JH, N
                DO KI=1,3
                        IK=3*(IH-1)+KI
                        JK=3*(JH-1)+KH
                        HESS(IK,JK)=(GRAD_PLUS(IH,KI)-GRAD_MIN(JH,KJ))/(2.0*DELTA)
                        HESS(JK,IK)=HESS(IK,JK)
                ENDDO
            ENDDO
            ENDDO
        ENDDO

! }}}

