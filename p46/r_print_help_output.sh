cat << EOF
OUTPUT OPTIONS
	varying output filenames:

	--append_job_number	append job number, only makes sense if the option -j is 
       					specified	
					set internal parameter append
	--n-append_job_number	don't append job number
	--append_job_name	append job name, specified by the switch --pbs_job_name 

	-od, --output_dir_name	specify the output directory name (default is out)
EOF
