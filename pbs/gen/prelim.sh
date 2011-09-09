
mkdir -p tmp
mkdir -p xyz

if [ ! $append_files ]; then 
	eo "-n" 
fi

start_date=` date_dmy_hms `
start_time=` date_in_secs `

eo "start date=$start_date" 	 
eo "start time=$start_time" 	 
eo


