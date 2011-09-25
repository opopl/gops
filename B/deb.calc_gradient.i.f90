       if (deb) then
!          write(fh,'(a)') '==============================='
          !write(fh,'(a,a20,i20)') 'calc_gradient ','Call ',NCALL
          !write(fh,'(a)') '==============================='
          !write(fh,'(A20,F20.5)') 'Total RMS =', RMS(1)
          !write(fh,'(A20,F20.5)') 'Non-bonded =', RMS(2)
          !write(fh,'(A20,F20.5)') 'bonded =', RMS(3)
          !write(fh,'(A20,F20.5)') 'bond angles =', RMS(4)
          !write(fh,'(A20,F20.5)') 'dihedral angles =', RMS(5)
          write(fh,101) NCALL, 'old_g',trim(ptype), (RMS(i),i=1,5)
          include "fmt.i.f90"
        endif

