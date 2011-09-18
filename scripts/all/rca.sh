#!/bin/bash

export shd="`dirname $(readlink -f $0)`"
export this_script=` basename $0 `
export gopsdir=$shd/../../

source $shd/rca.$1.i.sh
