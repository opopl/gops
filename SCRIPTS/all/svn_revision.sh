#!/bin/bash

# directory where this script resides
export shd="`dirname $(readlink -f $0)`"
# name of this script 
export this_script=` basename $0 `
export wg_dir="$shd/../../"
export inc_dir="$shd/../../INCLUDE"

dir="."
[ ! -z $1 ] && dir="$1" 
#svnversion $dir | sed 's/.*://' | sed 's/M//' > version.tmp

svn info $wg_dir >& n  
s=`cat n | awk '/Revision:/ { print $2 }'`
[ ! -z "$s" ] && echo "$s" > $inc_dir/SVNREV

if [ -f SVNREV ]; then
 	 cat $inc_dir/SVNREV
else
 	exit 1
fi
