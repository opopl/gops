cat << EOF
PBS 
	-q, --pbs_queue_name string	queue name
	
	--pbs_job_name 	 string		job name; default is the directory name 
					in which script is located

	-wt, --pbs_wall_time	string	specify PBS walltime
						  "day"	24:00:00
EOF

