module dv

implicit none 

contains 

subroutine display_version(file)

integer file
character(100) fflags(10)

write(file,10) "==========================================="
write(file,10)
write(file,10) "Gb - "
write(file,10) 
write(file,10) "Compilation time: Wed Sep 21 14:12:01 BST 2011"
write(file,10) "Compiled by op226@serenity on GNU/Linux x86_64"
write(file,10)
write(file,10) "Compiler executable:  gfortran"
! write(file,20) "Compiler flags:  ", "-ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic  -I."
write(file,10) "Command-line options passed to makefile:  "
write(file,10)
write(file,10) "==========================================="

10 format(a)
20 format(2a)

end subroutine

end module
