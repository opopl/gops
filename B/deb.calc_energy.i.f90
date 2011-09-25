       if (deb) then
          !write(fh,'(a)') '==============================='
          !write(fh,'(a,a20,i20)') 'calc_energy ','Call ',NCALL
          !write(fh,'(a)') '==============================='
          !write(fh,'(A20,F20.5)') 'Total BLN energy =', energy
          !write(fh,'(A20,F20.5)') 'Non-bonded =', e_nbond
          !write(fh,'(A20,F20.5)') 'bonded =', e_bond
          !write(fh,'(A20,F20.5)') 'bond angles =', e_bangle
          !write(fh,'(A20,F20.5)') 'dihedral angles =', e_tangle
          write(fh,101) NCALL, 'old_e',trim(ptype), (E(i),i=1,5)
          include "fmt.i.f90"
        endif

