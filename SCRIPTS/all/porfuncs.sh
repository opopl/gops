#!/bin/bash

# GPL License Info  {{{
#
#   Portability functions module generator
#   Copyright (C) 2003-2005 Semen A. Trygubenko and David J. Wales
#   This file is part of OPTIM.
#
#   OPTIM is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   OPTIM is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#}}}

# read in compiler info from command-line input {{{

compilers=( nag ifort pgi ifc g95 gfortran pathscale )

if [ $# -eq 1 ]; then 
  	compiler="$1"
  	case $compiler in	
    		nag|ifort|pgi|ifc|g95|gfortran|pathscale) ;;
		*) 
    			echo 'Unknown compiler!'
        		exit 1
		;;
  	esac
else
  	compiler="pgi"
fi

#}}}

echo "MODULE PORFUNCS" 
#{{{


case $compiler in 
 	nag)  
cat << EOF
use f90_unix, only: getarg
use f90_unix_proc, only: system, exit
EOF
;;
esac

cat << EOF
implicit none
contains
EOF

# FLUSH {{{
case $compiler in
  	nag) # {{{
cat << EOF
subroutine flush(UNIT,ISTAT) ! flushes the output buffer of logical unit UNIT which must be

use f90_unix, NAGflush => flush ! connected for formatted sequential output; ISTAT is ignored

implicit none

integer,intent(in) :: UNIT
integer,intent(out),optional :: ISTAT
call NAGflush(UNIT)

end subroutine flush
EOF
;;
#}}}
  	g95) #{{{
cat << EOF
subroutine flush(UNIT,ISTAT) ! flushes the output buffer of logical unit UNIT which must be

implicit none ! connected for formatted sequential output; ISTAT is ignored

integer,intent(in) :: UNIT
integer,intent(out),optional :: ISTAT

call flush(UNIT)

end subroutine flush
EOF
;;
#}}}
	gfortran) #{{{
cat << EOF
subroutine flush(UNIT,ISTAT) ! flushes the output buffer of logical unit UNIT which must be

implicit none ! connected for formatted sequential output; ISTAT is ignored

integer,intent(in) :: UNIT
integer,intent(out),optional :: ISTAT
!              call flush(UNIT)

end subroutine flush
EOF
;;
# }}}
esac
 
#}}}
# GETARG (added by op226) {{{

echo "          subroutine getarg_subr(position,value) ! wraps getarg function so it can be use-associated"

#nag{{{
case $compiler in
  	nag) echo "               use f90_unix, only: getarg" ;;
esac
#}}}

cat << EOF

implicit none

integer,intent(in) :: position
character(len=*),intent(out) :: value

call getarg(position,value)

end subroutine getarg_subr
EOF


#}}}
# IARGC {{{
cat << EOF

subroutine iargc_subr(n) ! wraps iargc function so it can be use-associated

EOF

case $compiler in
  	nag) echo "use f90_unix, only: iargc" ;;
esac

cat << EOF

implicit none
integer,intent(out) :: n

EOF


if [[ ! $compiler == "nag" ]]; then
echo "integer iargc"
fi

cat << EOF
n = iargc()
end subroutine iargc_subr

EOF
 
#}}}
# FORK {{{
echo "          subroutine fork_subr(pid)" ! returns zero in the child process, PID of child in parent process

case $compiler in
  	nag) echo "               use f90_unix_proc, only: fork" ;;
esac

echo "               implicit none"
echo "               integer, intent(inout) :: pid"
echo " "
if [ $compiler == "pgi" ]; then
     echo "               integer fork"
fi

echo " "

case $compiler in
  	pgi) 	echo "               pid=fork()" ;;
      ifort) 	echo "               integer ierror"
     		echo "               call pxffork(pid,ierror)"
     	;;
	nag) 	echo "               call fork(pid)";;
	pathscale) 
		echo "               integer fork"
     		echo "               pid=fork()"
		;;
	gfortran);;
esac

echo "          end subroutine fork_subr"
echo " "
#}}}
# SYSTEM {{{
echo "          subroutine system_subr(JobString,ExitStatus)"

if [ $compiler == 'nag' ]; then
	  echo "               use f90_unix_proc, only: system"
fi

cat << EOF
               implicit none
 
               character(len=*),intent(in) :: JobString
               integer,intent(out) :: ExitStatus
EOF
 
case $compiler in 
	ifort) # {{{
		  echo "               integer shiftr,system"
		  echo " "
		  echo "               ExitStatus=system(JobString)"
		  ;;
		#}}}
	pgi) # {{{
		  echo "               integer system"
		  echo " "
		  echo "               ExitStatus=system(JobString)"
		  echo "               ExitStatus=ishft(ExitStatus,-8)"
		  ;;
		  # }}}
	nag) # {{{
		  echo "               call system(JobString,ExitStatus)"
		  echo "               ExitStatus=ishft(ExitStatus,-8)"
		  ;;
		#}}}
	pathscale) # {{{
		  echo "               integer system"
		  echo "               ExitStatus=system(JobString)"
		  ;;
		#}}}
esac

cat << EOF
          end subroutine system_subr
 
EOF

#}}}
# WAIT {{{
echo "          subroutine wait_subr(pid,ExitStatus)"
if [[ $compiler == 'nag' ]]; then
  echo "               use f90_unix_proc, only: wait"
fi

cat << EOF
               implicit none
 
               integer,intent(inout) :: pid,ExitStatus
EOF
 
case $compiler in
ifort) # {{{
   echo "               integer shiftr,ierror"
   echo " "
   echo "               call pxfwait(ExitStatus,pid,ierror)"
   ;;
 #}}}
pgi) # {{{
   echo "               integer wait"
   echo " "
   echo "               pid=wait(ExitStatus)"
   echo "               ExitStatus=ishft(ExitStatus,-8)"
   ;;
 #}}}
nag) #{{{
  echo "               call wait(ExitStatus,pid)"
  echo "               ExitStatus=ishft(ExitStatus,-8)"
  ;;
  # }}}
pathscale) # {{{
  echo "               INTEGER WAIT"
  echo "               pid=wait(ExitStatus)"
  ;;
  # }}}
  gfortran) ;;
esac

echo "end subroutine wait_subr"
 
#}}}
# GETPID {{{
echo "          subroutine getpid_subr(pid)"
if [[ $compiler == 'nag' ]]; then
  echo "               use f90_unix, only: getpid"
fi

cat << EOF
               implicit none
 
               integer,intent(out) :: pid
EOF

case $compiler in
  	ifort|pgi|ifc|pathscale)
  echo "               integer getpid"
  echo " "
  ;;
esac

cat << EOF
               pid=getpid()
          end subroutine getpid_subr
EOF
#}}}
#}}}
echo "END MODULE PORFUNCS"
