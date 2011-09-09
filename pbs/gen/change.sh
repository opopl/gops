#!/bin/bash

# input file with keywords

f=data.G46

while [ ! -z "$1" ]
	do 
		case "$1" in
	"f" | "force") cat $f | sed "/PULL/ s/[0-9\.DE-]*$/$2/" >& n; mv n $f ;;
	"t" | "temperature") 
		awk '/TEMPERATURE/  { sub("^[0-9.]*", par, $2);  }{ print $0 }' par=$2 $f > n ; mv n $f ;;
	"nsteps")      awk '/STEPS/  { sub("^[0-9ED.]*", par, $2);  }{ print $0 }' par=$2 $f > n; mv n $f ;;
	"target")      awk '/TARGET/  { sub("^[0-9ED.-]*", par, $2);  }{ print $0 }' par=$2 $f > n; mv n $f ;;
	"pull_start")  awk '/PULL/  { sub("^[0-9.]*", par, $2);  }{ print $0 }' par=$2 $f > n; mv n $f ;;
	"pull_end")    awk '/PULL/  { sub("^[0-9.]*", par, $3);  }{ print $0 }' par=$2 $f > n; mv n $f ;;
	"radius")      awk '/RADIUS/  { sub("^[0-9.]*", par, $2);  }{ print $0 }' par=$2 $f > n; mv n $f ;;
	"system") ;;
		esac
	shift
done
