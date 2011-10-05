
      SUBROUTINE FINALIO
!op226> Declarations {{{ 
      USE COMMONS
      USE V
      USE F

      IMPLICIT NONE

      INTEGER J1, J2, J3
      character(len=30) t
      ! }}}
! subroutine body {{{

      OPEN(LE_FH,FILE=LE_FILE,STATUS='UNKNOWN')
      OPEN(EA_FH,FILE=EA_FILE,STATUS='UNKNOWN')
      OPEN(FOFH,FILE=FO_FILE,STATUS='UNKNOWN')

      CALL GETTIME(T)
      T="# Time: "//adjustl(T)

      call ed(fofh)
      write(fofh,'(a)') 'Final Output'
      write(fofh,'(a)') T 
      call ed(fofh)
      write(fofh,'(a)') 'Energies'
      write(fofh,'(a)') ''
      write(fofh,'(a)') 'Global minimum: ' 
      write(fofh,'(a)') ''
      write(fofh,'(2a20,1x,e20.5)') 'Total energy: ', ' E_total=', eamin(1,1)
      write(fofh,'(2a20,1x,e20.5)') 'Non-bonded interaction energy: '  ,'E_nbond= ', eamin(1,2)
      write(fofh,'(2a20,1x,e20.5)') 'Bonded interaction energy: ','E_bond= ', eamin(1,3)
      write(fofh,'(2a20,1x,e20.5)') 'Bond angle interaction energy: ','E_bangle= ', eamin(1,4)
      write(fofh,'(2a20,1x,e20.5)') 'Torsional angle interaction energy: ','E_tangle= ', eamin(1,5)
      if (pullt) then 
         write(fofh,'(2a20,1x,e20.5)') 'Applied force energy: ','E_dfz= ', eamin(1,6)
      endif
      write(fofh,'(2a20,1x,e20.5)') 'Radius of gyration: ',' R_gyr=', rgmin(1)
      write(fofh,'(a)') ''
      write(fofh,'(i5,2x,a)') NSAVE,'Saved lowest energies:' 
      write(fofh,'(a)') ''
      write(fofh,'(e20.5)') eamin(1:nsave,1) 
      write(fofh,'(a)') ''
      call ed(fofh)

      WRITE(EA_FH,'(A1,A)') "# Command-line: ",CMDLINE
      WRITE(EA_FH,'(A)') T
      WRITE(EA_FH,'(E10.5,6E20.5)') PFORCE,EAMIN(1,1:6)

      DO J1=1,NSAVE
      ! {{{

         WRITE(LE_FH,'(I10)') NATOMS
         WRITE(LE_FH,10) J1,QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)
            DO J2=1,NATOMS
               WRITE(LE_FH,'(1A1,1X,3F20.10)') BEADLETTER(J2),(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
       ENDDO

      CLOSE(LE_FH)
      CLOSE(EA_FH)
      CLOSE(FOFH)

      ! }}}
      RETURN
      END
     
