
! Specify amino acid types by filling in array NTYPE(:)
! 1 => Hydrophobic (B);  2 => Hydrophilic (L); 3 => Neutral (N)
        NTYPE(1:9) = 1
        NTYPE(10:12) = 3
        NTYPE(13:19:2) = 2
        NTYPE(14:20:2) = 1
        NTYPE(21:23) = 3
        NTYPE(24:32) = 1
        NTYPE(33:35) = 3
        NTYPE(36:46:2) = 2
        NTYPE(37:45:2) = 1

