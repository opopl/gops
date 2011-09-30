#!/usr/bin/env bash

# set_base_vars set_base_dirs display_help {{{

set_base_dirs(){
# {{{
export packdir=$HOME/arch/packed
export unpackdir=$HOME/arch/unpacked
export gopsdir="$shd/../../"
# directory where this script resides
export shd="`dirname $(readlink -f $0)`"
# }}}
}

set_base_vars(){
# {{{
s_purpose="generate compiler-specific makefile-include files"
s_project="gops/scripts/all"

# name of this script 
export this_script=` basename $0 `

vim_opts="-n -p"
v="vim $vim_opts"
repos=( "config" "scripts" "templates" "vrt" "install" "doc-coms" "doc-cit" )

# }}}
}

set_base_vars

display_help(){
# {{{
# general beginning part {{{
cat << EOF
=============================================
SCRIPT NAME: $this_script 
PROJECT: $s_project
PURPOSE: $s_purpose
USAGE: $this_script [ COMPILER ] [ OPTIONS ] 

	OPTIONS:

	============
	General
	============

			display the help message

	vm		v(iew) m(yself), i.e., edit this script
	-g		gvim invocation
	
	============
EOF
# }}}
# actual interesting part {{{
cat << EOF

	COMPILERS:
		gf (gfortran)
		nag 
		ifort
		pgi
	OPTS:
		debug 
		opt

EOF
# }}}
# final general part {{{
cat << EOF
REMARKS:
AUTHOR: O. Poplavskyy
SCRIPT LOCATION:
    $shd/$this_script
=============================================
EOF
# }}}
# }}}
}

[ -z "$*" ] && ( display_help; exit 0 )

#}}}

main(){
# {{{

while [ ! -z $1 ]; do   
	case "$1" in
	  pgi|gf|gfortran|nag|ifort|gf) #{{{
#define FC, intro part {{{
compiler=$1
case "$compiler" in
  	pgi) FC=pgf90 ;;
  	nag) FC=nagfor ;;
  	ifort) FC=ifort ;;
  	gf) FC=gfortran ;;
esac

cat << EOF
#
# comp.mk. Compiler include file
# Generated: `date` 
# Compiler: $compiler
#
export FC=$FC
EOF
cat << 'EOF'
# path for the objects
export OBJSPATH=$(ROOTPATH)/obj/$(PROGNAME)/fc_$(FC)/
# modules path for specific compiler
export MODPATH=$(MODGENPATH)/fc_$(FC)/
# library path for *.a compiled with $(FC)
export LIBAPATH=$(LIBPATH)/fc_$(FC)/
EOF
#}}}
case "$1" in
  # pgi {{{
  pgi) 
cat << 'EOF'
MODFLAG:= -I. -module $(MODPATH)
FFLAGS_o := -fast -Mipa=fast,inline -Msmartalloc
FFLAGS_g := -Mextend -O0 -Mnoframe -g -traceback
LDFLAGS= -L.
SWITCH=pgi
EOF
;;
#}}}
# nag {{{
nag)
cat << 'EOF'
MODFLAG= -mdir $(MODPATH)
FFLAGS_g = -132 -maxcontin=3000 -C=all -mtrace=all -gline -kind=byte 
FFLAGS_g = -132 -g90 
LDFLAGS= -L.
SWITCH=nag
EOF
;;
#}}}
# ifort {{{
ifort)
cat << 'EOF'
FFLAGS_g= -g -debug all -check all -implicitnone -warn unused -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds
FFLAGS_o= -132 -O4
SWITCH=ifort
LDFLAGS= -L.
MODFLAG=-I. -module $(MODPATH)
EOF
;;
#}}}
# gfortran {{{
gf|gfortran) 
cat << 'EOF'
MODFLAG=-I. -J$(MODPATH)
SWITCH=gfortran
LDFLAGS =-L.
FFLAGS_g= -ffixed-line-length-none -g -fbounds-check -Wuninitialized -O -ftrapv 
#FFLAGS_g += -fimplicit-none -fno-automatic 
FFLAGS_o= -ffixed-line-length-none -O3 -ftree-vectorize
EOF
;;
#}}}
esac
cat << 'EOF'
FFLAGS = $(FFLAGS_g) 
EOF
	;;
#}}}
	debug|opt)
		case $1 in
# debug {{{
debug)
cat << 'EOF'
FFLAGS = $(FFLAGS_g) 
EOF
;;
#}}}
# opt {{{
opt)
cat << 'EOF'
FFLAGS = $(FFLAGS_o) 
EOF
;;
#}}}
		esac	
		echo 'FFLAGS+=$(MODFLAG)'
	;;
	esac
	shift
done 

# }}}
}

# main part 
# {{{

script_opts=( $* )
set_base_dirs

while [ ! -z "$1" ]; do
  	case "$1" in
		  #{{{
	  	vm) $v $0 $hm/scripts/f.sh; exit ;;
		h) display_help $*; exit ;;
		-g) v="$v -g" ;;
	  	*) main $* && exit 0 ;;
	esac
  	shift
        #}}}
done

# }}}


