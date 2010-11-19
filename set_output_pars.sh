
if [ -f job_number ]; then
	let job_number=` cat job_number ` 
else
  	job_number=""
fi

export output_dir=$data_dir/$base_output_dir
export main_log_file=$output_dir/r.log$job_number

ef_header_line="$datahdr $start_date; \
	N=`printvar 'nsteps' $nsteps` ; \
	F=`printvar 'force' $force`; T=$temp; nrg=$nrg; R=$radius; \
	j=$pbs_job_id; s=$sys"

# log file r.log is not included

base_output_files=( gm.xyz ef.dat pdff.dat gmin_out hits pdf energy.dat markov.dat tl.dat )

let i=0

for f in ${base_output_files[@]}
	do
	  	output_files[$i]=$output_dir/$f$job_number
		case "$f" in
		  	"ef.dat") output_file_lines[$i]="$ef_header_line" ;;
		esac
		let i=$(($i+1))
done

export plot_dir="$data_dir/plots"
