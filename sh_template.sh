#!/bin/bash
  
source sh_header.sh

if [ -z "$*" ]; then 
  	source "$this_script"_print_help.sh
else

while [ ! -z "$1" ] 
	do
	  source "$this_script"_read_cmd.sh
	  shift
done 

source "$this_script"_ex.sh

fi

