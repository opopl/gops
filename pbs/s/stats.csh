#!/bin/csh

set radius=$1
set natoms=46
# set exec=~/svn/GMIN/bin/GMIN.4.0
set exec=GMIN

set directory=$2
set seed=$3
set nrg=$4

mkdir -p $directory
cd $directory
cp ../data.G46 data
rm hits >& /dev/null
set count=1

echo
echo Running GMIN for P46 from $nrg random starting geometries
echo Random number seeds start from $seed
echo

while ($count <= $nrg)
   echo $natoms $radius -$seed > randata
   rancoords # >& /dev/null
   cp newcoords coords
   echo count $count
   $exec >& $directory.output.$natoms.$count
   echo `grep hit GMIN_out | head -1 | sed -e 's/[a-zA-Z]//g' -e 's/[a-zA-Z]//g' -e 's/\.//' -e 's/>//'` \
        `grep time= GMIN_out | tail -1 | sed 's/.*time=//'` >> hits
   cp coords coords.$count
   @ count +=1
   @ seed +=1 
end

gminconv2 < hits > temp ; head -1 temp > pdf
echo
echo Mean and standard deviation for global minimum first encounter time
echo
cat pdf
