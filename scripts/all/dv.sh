#!/bin/bash

export shd="`dirname $(readlink -f $0)`"
export this_script=` basename $0 `

w0="write(file,10)"

vars=( "prog" "fflags" "fc_exec" "make_opts" )

# get parameters, e.g, compiler etc., from command line 
# {{{
while [ ! -z "$*"  ]; do
  	var="$1"
	for known_var in ${vars[@]}
		do
			case "$var" in
			  	"$known_var") export $var="$2" ;;
			esac
        done
	shift
done

case "$prog" in
  	G|gmi|gmin) prog_full="A program for finding global minima" ;;
	opti|optim) prog_full="A program for optimizing geometries and calculating reaction pathways" ;;
	pths|pathsample) prog_full="A driver for optim to create stationary point databases and perform kinetic analysis"
	;;
esac
# }}}

# send to output: subroutine display_version(file)
# {{{
cat << EOF      
module dv

implicit none 

contains 

subroutine display_version(file)

integer file
character(100) fflags(10)

$w0 "==========================================="
$w0
$w0 "$prog - $prog_full"
$w0 
$w0 "Compilation time: ` date `"
$w0 "Compiled by $USER@$HOSTNAME on ` uname -o` ` uname -m`"
$w0
$w0 "Compiler executable:  $fc_exec"
! write(file,20) "Compiler flags:  ", "$fflags"
$w0 "Command-line options passed to makefile: $make_opts "
$w0
$w0 "==========================================="

10 format(a)
20 format(2a)

end subroutine

end module
EOF
# }}}
