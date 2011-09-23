
      include "bln.ntype.inc.f90"         ! specify bead types 

        SELECTCASE(PTYPE)
                CASE("GO")                ! Go-like model
      ! ==================================================
      include "bln.go.connect.inc.f90"    ! specify native contacts 
      include "bln.go.ab.inc.f90"         ! parameters for the non-bonded LJ interaction 
      include "bln.go.cd.inc.f90"         ! parameters for the torsion angle interaction
      ! ==================================================
                CASE("WT")                ! Wild-type (WT) original frustrated BLN model 
      ! ==================================================
      include "bln.wt.ab.inc.f90"         ! L-J interaction between non-bonded particles
      include "bln.go.cd.inc.f90"         ! Torsion angle potential
      ! ==================================================
        ENDSELECT

