
        do i=2,N-1
                AN=ANG(I,1)
                GBA_I(i,1:3)=-(F(i,1)/(B(i-1)*sin(AN)))*(EB(i,1:3)+EB(i-1,1:3)*cos(AN))
                GBA_K(i,1:3)=(F(i,1)/(B(i)*sin(AN)))*(EB(i,1:3)*cos(AN)+EB(i-1,1:3))
                GBA_J(i,1:3)=-GBA_I(i,1:3)-GBA_K(i,1:3)
        enddo

        GBA(1,1:3)=GBA_I(2,1:3)
        GBA(2,1:3)=GBA_I(3,1:3)+GBA_J(2,1:3)

        do i=3,N-2
                GBA(i,1:3)=GBA_I(i+1,1:3)+GBA_J(i,1:3)+GBA_K(i-1,1:3)
        enddo 

        GBA(N-1,1:3)=GBA_J(N-1,1:3)+GBA_K(N-2,1:3)
        GBA(N,1:3)=GBA_K(N-1,1:3)
