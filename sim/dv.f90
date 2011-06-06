module dv

implicit none 

contains 

subroutine display_version(file)

integer file

write(file,10) "==========================================="
write(file,10)
write(file,10) "GMIN - A program for finding global minima"
write(file,10) 
write(file,10) "Compilation time: Mon Jun  6 19:46:04 BST 2011"
write(file,10) "Compiled by gates@Gates-PC5 on GNU/Linux i686"
write(file,10)
write(file,10) "Compiler executable:  gfortran"
!write(file,10) "Compiler flags: -ffree-form -I."
write(file,10) "Command-line options passed to makefile:  "
write(file,10)
write(file,10) "==========================================="

10 format(a)

end subroutine

end module
