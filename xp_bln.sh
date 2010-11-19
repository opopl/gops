#!/bin/bash

# bln.sh: this will give array xyz_headers 
# 		 plus insert correct atom labels into the xyz file
#		 plus it removes blank lines 

#
# the original BLN sequence is
# is B9-N3-(LB)4-N3-B9-N3-(LB)5-L 
#


# remove blank lines

cat $full_file | sed '/^$/d' > n; mv n $full_file

subs(){

fname=$4
gawk -f $shd/"$this_script"_sub_symbol.awk -v start=$1 -v end=$2 -v symbol="$3" $fname >> $fname.1

}

f="$full_file.n"
cp $full_file $f

file_size=$(stat -c%s "$f")

fnew=` mktemp `

let i=0
while [ $file_size -gt 0 ]; do

	xyz=` mktemp `
	head=` mktemp `
	n=` mktemp ` 
  	
	cat $f | sed -n '1,2p' >& $head
	xyz_header=` cat $f | sed -n "2p" `

  	cat $f | sed -n '3,48p' >& $xyz
	
	cat $f | sed '1,48d' >& $n ; mv $n $f
	
	rm -f $xyz.1
	
	# B-9
	subs 1 9 "H" $xyz
	# N-3
	subs 10 12 "N" $xyz
	# (LB)-4
	for (( j=0;j<4;j++)); do	
		let j1=$((13+2*$j))
	  	let j2=$((14+2*$j))
		subs $j1 $j1 "C" $xyz
		subs $j2 $j2 "H" $xyz
      	done
	# N-3
	subs 21 23 "N" $xyz
	# B-9
	subs 24 32 "H" $xyz
	# N-3
	subs 33 35 "N" $xyz
	# (LB)-5
	for (( j=0;j<5;j++)); do	
		let j1=$((36+2*$j))
	  	let j2=$((37+2*$j))
		subs $j1 $j1 "C" $xyz
		subs $j2 $j2 "H" $xyz
      	done
	# L
	subs 46 46 "C" $xyz

	cat $head >> $fnew
	cat $xyz.1 >> $fnew

	file_size=$(stat -c%s "$f")
	let i=$(($i+1))
	
	xyz_headers[$i]=$xyz_header
	forces[$i]=` echo $xyz_header | awk '{ print $4 }' | sed 's/;$// ' ` 

	rm -r $head $xyz $xyz.1
done

let num_frames=$i

cp $fnew $full_file
