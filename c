#!/bin/bash

dir_name=` basename $PWD `
this_script=` basename $0 ` 

ods=( out dn )

if [ -z "$1" ]; then  
cat << EOF
============================================
SCRIPT NAME:	$this_script
PURPOSE:	clean files
USAGE: $this_script OPTS
OPTS:
		display this help message
	a	all
	p	plots
	e	pbs error files
	q	q.log only
	dn	rm -rf dn-*
	o	output dirs: out* dn*
	* 	rm \$1*
============================================
EOF
	else
	  	while [ ! -z "$1" ]; do  
		  	case "$1" in
				e) rm -f $dir_name.o*  ;;
				a) $this_script e q p dn o t ef ;;
				p) rm -f plots/*.* ;;
				q) echo "" >& q.log ;;
				dn) rm -rf dn-* ;;
				o) for od in ${ods[@]}
			       		do	
						 rm -rf $od*
			      	done
				;;
				t) rm -rf t.dat ;;
				ef) cd ef ; rm -f *.pdf *.ps ; cd - ;;
				*) rm -rf $1* ;;
			esac
			shift
	done
fi

#vim:set ft=sh
