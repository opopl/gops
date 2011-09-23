module dv

implicit none 

contains 

subroutine display_version(file)

integer file
character(100) fflags(10)

write(file,10) "==========================================="
write(file,10)
write(file,10) "B - "
write(file,10) 
write(file,10) "Compilation time: Fri Sep 23 23:48:49 BST 2011"
write(file,10) "Compiled by op226@mek-quake on GNU/Linux x86_64"
write(file,10)
write(file,10) "Compiler executable:  pgf90"
! write(file,20) "Compiler flags:  ", "-Mextend -O0 -Mnoframe  -module /home/op226/gops/B/..//mod/B//fc_pgf90/ "
write(file,10) "Command-line options passed to makefile:  "
write(file,10)
write(file,10) "==========================================="

10 format(a)
20 format(2a)

end subroutine

end module
