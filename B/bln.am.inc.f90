
        ALLOCATE(R(N,3),GRAD(N,3))
        ALLOCATE(CONNECT(N,N))

        ALLOCATE(AB(N,N,2))
        ALLOCATE(CD(N,2))

        allocate(DR(N,N,3))
        allocate(LEN_DR(N,N))

        allocate(ANG(N,2))
        allocate(F(-1:N+1,2))
        allocate(F(N,2))
