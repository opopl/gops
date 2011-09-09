
#cat $qn | sed "s/^\#PBS -N [ .]*/\#PBS -N $job_name/ " > $qnm
awk '/\#PBS -q/ { sub("^[a-z0-9]*", name , $3) } { print $0 }' name=$queue_name $qn > n; mv n $qnm
