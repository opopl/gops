        ! a=1, 2, 3 and N-2, N-1, N are special {{{

        ! a=1
        GTA(1,1:3)=-PP(1,1:3)*F(2,2)*B(2)
        ! a=2
        GTA(2,1:3)=PP(2,1:3)*(-FB(3,1)+DPD(2)*FB(2,2)) &
                  +PP(1,1:3)*( FB(2,1)-DPD(1)*FB(2,2))
        ! a=3
        GTA(3,1:3)=PP(3,1:3)*(-FB(4,1)+DPD(3)*FB(3,2)) &
                  +PP(2,1:3)*( FB(3,1)-FB(2,1)-DPD(2)*(FB(3,2)+FB(2,2))) &
                  +PP(1,1:3)*DPD(1)*FB(2,2)
        ! }}}
        
        ! a=4,...,N-3
        do i=4,N-3
            GTA(i,1:3)=PP(i,1:3)*(-FB(i+1,1)+DPD(i)*FB(i,2)) &
                  +PP(i-1,1:3)*( FB(i,1)-FB(i-1,1)-DPD(i-1)*(FB(i,2)+FB(i-1,2))) &
                  +PP(i-2,1:3)*( DPD(i-2)*FB(i-1,2) + FB(i-2,1) )
        enddo

        ! a=N-2,N-1,N {{{

        ! a=N-2
        GTA(N-2,1:3)=PP(N-2,1:3)*DPD(N-2)*FB(N-2,2) &
                    +PP(N-3,1:3)*( FB(N-2,1)-FB(N-3,1)-DPD(N-3)*(FB(N-2,2)+FB(N-3,2))) &
                    +PP(N-4,1:3)*( DPD(N-4)*FB(N-3,2) + FB(N-4,1) )
        ! a=N-1
        GTA(N-1,1:3)= &
                  PP(N-2,1:3)*( FB(N-2,1)-DPD(N-2)*FB(N-2,2)) &
                  +PP(N-3,1:3)*( DPD(N-3)*FB(N-2,2) + FB(N-3,1) )
        ! a=N
        GTA(N,1:3)=PP(N-2,1:3)*FB(N-2,1)
        ! }}}
