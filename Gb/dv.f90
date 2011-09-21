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
write(file,10) "Compilation time: Wed Sep 21 22:01:08 BST 2011"
write(file,10) "Compiled by op226@mek-quake on GNU/Linux x86_64"
write(file,10)
write(file,10) "Compiler executable:  pgf90"
! write(file,20) "Compiler flags:  ", "-Mextend -O0 -Mnoframe "
write(file,10) "Command-line options passed to makefile:  "
write(file,10)
write(file,10) "==========================================="

10 format(a)
20 format(2a)

end subroutine

end module
