
export pref="$pref_default"

#if [ $pref_pbs_job_id -eq 1 ]; then
          #export pref="$pref"."$pbs_job_id" 
#fi 

#if [ $pref_pbs_job_name -eq 1 ]; then
          #export pref="$pref"."$job_name" 
#fi 

#if [ $pref_job_number -eq 1 ]; then
          #export pref="$pref"."$job_number" 
#fi 

#source f.sh
