MODULE PORFUNCS
! Generated on: Wed Sep 28 01:00:21 BST 2011 
! Compiler: ifort
implicit none
contains
          subroutine getarg_subr(position,value) ! wraps getarg function so it can be use-associated

implicit none

integer,intent(in) :: position
character(len=*),intent(out) :: value

call getarg(position,value)

end subroutine getarg_subr

subroutine iargc_subr(n) ! wraps iargc function so it can be use-associated


implicit none
integer,intent(out) :: n

integer iargc
n = iargc()
end subroutine iargc_subr

          subroutine fork_subr(pid) ! returns zero in the child process, PID of child in parent process
               implicit none
               integer, intent(inout) :: pid
 
 
               integer ierror
               call pxffork(pid,ierror)
          end subroutine fork_subr
 
          subroutine system_subr(JobString,ExitStatus)
               implicit none
 
               character(len=*),intent(in) :: JobString
               integer,intent(out) :: ExitStatus
               integer shiftr,system
 
               ExitStatus=system(JobString)
          end subroutine system_subr
 
          subroutine wait_subr(pid,ExitStatus)
               implicit none
 
               integer,intent(inout) :: pid,ExitStatus
               integer shiftr,ierror
 
               call pxfwait(ExitStatus,pid,ierror)
end subroutine wait_subr
          subroutine getpid_subr(pid)
               implicit none
 
               integer,intent(out) :: pid
               integer getpid
 
               pid=getpid()
          end subroutine getpid_subr
END MODULE PORFUNCS
