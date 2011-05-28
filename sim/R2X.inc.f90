      DO K=1,3
        IXMIN=1+(K-1)*NA
        IXMAX=K*NA
        X(IXMIN:IXMAX)=R(1:NA,K)
      ENDDO

