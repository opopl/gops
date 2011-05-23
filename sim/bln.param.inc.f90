
       include bln.ntype.inc.f90        ! specify bead types 

        SELECTCASE(PTYPE)
                CASE("GO")              ! Go-like model
      include bln.go.connect.inc.f90    
      include bln.go.ab.inc.f90         
      include bln.go.cd.inc.f90
                CASE("WT")              
                        include bln.wt.ab.inc.f90       ! L-J interaction between non-bonded particles.
                        include bln.go.cd.inc.f90       ! Dihedral angle potential
        ENDSELECT

