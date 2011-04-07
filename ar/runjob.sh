#!/bin/sh
#
#-------- SGE parameters ----------
#$ -S /bin/bash
#$ -N test -o test.out -e qsub.out
#$ -V -cwd -j y -m n -notify
#$ -M jonas@chem.helsinki.fi
##$ -l chem=true
#$ -pe mpi 3
#-------- end SGE parameters ----------
#
#-------- PBS parameters ----------
#PBS -o run.out -e qsub.out
#PBS -q itanium
##PBS -lnodes=4:ppn=4
#PBS -lnodes=1
#PBS -lwalltime=00:10:00
#PBS -lpmem=500MB
##PBS -m abe
#-------- end PBS parameters ----------

##
##   Writen by Jonas Juselius <jonas@iki.fi> 
##

###------ JOB SPECIFIC DEFNITIONS ------###

# possible jobs types are: mpich, lam, scali or serial
jobtype="serial" 
stop_on_crash=yes

# setup environment
#. $HOME/bin/turbomole.sh 

# a job id is automatically added to the workdir
workdir=/work/$USER 
bindir="$TURBODIR/bin/ia64-unknown-linux-gnu"
scriptdir="$HOME/scripts"

# name of program
executable="ridft"
# inputfiles are copied to the nodes, but not copied back
inputfiles="" 
# inoutfiles are copied to the nodes and back when the job is complete
inoutfiles="*"

outputdir="out"
copy_executable_to_nodes="no" # run over NFS or not

# command to run program
jobcmd="$executable"

# script(s) to (selectively) setup or clean unwanted files, and other stuff
prejob="fix_this_and_that.sh" 
postjob="cleanup.sh"

## setup copyall feature. copyall copies everything from all nodes back
## to the submit directory
do_copyall=no
copyall_tag_files=yes
copyall_max_size=10000   # size in kB

### Additional environment variables ###

###------ NO MODIFICATIONS NEEDED BEYOND THIS LINE  ------###

## begin function definitions 
setup_work_env()
{
	if test $PBS_O_WORKDIR; then
		submitdir=$PBS_O_WORKDIR
		jobid=$PBS_JOBID
		nodefile=$PBS_NODEFILE
	elif test $SGE_O_WORKDIR; then
		submitdir=$SGE_O_WORKDIR
		jobid=$JOB_ID
		make_sge_nodefile
	else
		submitdir=`pwd`
		jobid=$$
	fi
	if test "$copy_executable_to_nodes" = "yes"; then
		runbindir="$workdir"
	else 
		runbindir="$bindir"
	fi
	if test "x$outputdir" = "x"; then
		outputdir="$submitdir"
	else
		cwd=`pwd`
		cd $submitdir
		test ! -d "$outputdir" && mkdir $outputdir
		cd $outputdir
		outputdir=`pwd`
		cd $cwd
	fi
	workdir="$workdir/$jobid"
	error_flag=0
}

setup_workcopyrm()
{
	local scripts=""
	for i in $prejob $postjob; do
		if test -x $scriptdir/$i; then
			scripts="$scripts $scriptdir/$i"
		fi
	done
	if test "$copy_executable_to_nodes" = "yes"; then
		workcopy="$bindir/$executable $inputfiles $inoutfiles $scripts"
		workrm="$executable $inputfiles $scripts"
	else
		workcopy="$inputfiles $inoutfiles $scripts"
		workrm="$inputfiles $scripts"
	fi
}

make_sge_nodefile() 
{
    nodefile=/tmp/mpinodes-$USER.$HOSTNAME
	cat /dev/null > $nodefile
	cat $PE_HOSTFILE | while read line; do
		host=`echo $line | cut -f1 -d" "| cut -f1 -d"."`
		nslots=`echo $line | cut -f2 -d" "`
		echo "${host} cpu=${nslots}" >> $nodefile
	done
	cat $nodefile
	has_nodefile=1
}

start_lam()
{
	lamboot -ssi boot rsh -ssi rsh_agent "ssh -q -x" $nodefile
}

rm_nodefile()
{
	if test $has_nodefile; then
		test -f $nodefile && rm -f $nodefile
	fi
}

runcmd()
{
	case $jobtype in
		lam) lamexec N $* ;;
		mpich) mpiexec -comm none -pernode $* ;;
		scali) mpiexec -comm none -pernode $* ;;
		serial) $* ;;
	esac
}

dump_jobinfo()
{
	touch $1
	for i in `uniq $nodefile`; do
		echo $i:$workdir >>$1
	done
}

execjob()
{
	case $jobtype in
		lam) mpirun C $* ;;
		mpich) mpiexec $* ;;
		scali) mpimon $* -- `cat $nodefile` ;;
		serial) $* ;;
	esac

	if test "$?" != "0"; then
		test $stop_on_crash = yes && error_flag=1
	fi
}

make_crash_script()
{
local crash_script="$submitdir/cleanup_crash.$jobid"

echo "#!/bin/sh" > $crash_script
for i in `uniq $nodefile`; do
	echo "echo \"ssh $i rm -rf $workdir\"" >>$crash_script
	echo "ssh $i rm -rf $workdir" >>$crash_script
done
chmod 700 $crash_script
}

# generate a script in $workdir to facilitate copy with globbing 
# it's a bit primitive, and does not handle directories  very well
setup_copyall_script()
{
local tmp_script="$submitdir/copyall.sh.$$"
cat << EOF >$tmp_script
#!/bin/sh
rm -f $workdir/copyall.sh
ext=""
test $copyall_tag_files = yes && ext=".$HOSTNAME"
cd \$1
for i in *; do
	test \$i = "*" && exit 0
	test -z \$i && exit 0
	if test -d \$i; then 
		sz=\`du -ks \$i | sed -e 's/ *\([0-9]*\)\(.*\)/\1/'\`
		if (( $copyall_max_size < 0 )); then
			cp -r \$i \${2}/\${i}\${ext}
		elif (( \$sz < $copyall_max_size )); then 
			cp -r \$i \${2}/\${i}\${ext}
		fi
	else
		sz=\`du -k \$i | sed -e 's/ *\([0-9]*\)\(.*\)/\1/'\`
		if (( $copyall_max_size < 0 )); then
			cp \$i \${2}/\${i}\${ext}
		elif (( \$sz < $copyall_max_size ));then 
			cp \$i \${2}/\${i}\${ext}
		fi
	fi
done
EOF
chmod 700 $tmp_script
runcmd cp $tmp_script $workdir/copyall.sh
rm -f $tmp_script
}

finalize_job()
{
	test "$jobtype" = "lam" && lamhalt
	rm_nodefile
	:
}

## end functions

##
## Run the job
##
setup_work_env
setup_workcopyrm
cd $submitdir

## start LAM-MPI
test "$jobtype" = "lam" && start_lam

## create local workdir on every node
runcmd mkdir -p $workdir

for i in $workcopy; do
	[ -f $i ] && runcmd cp -r $i $workdir
done

# run any pre job scritp(s)
for i in $prejob; do
	[ -f $i ] && runcmd $workdir/$i 
done

cd $workdir

# Do the thing you do, baby!
if test -x "$runbindir/$executable"; then 
	execjob $runbindir/$jobcmd
else
	echo "Error: Could not execute $runbindir/$jobcmd"
fi

# if the job crashed, dump job info and exit for post mortem
if test $error_flag != 0; then 
	make_crash_script
	dump_jobinfo $submitdir/crash.$jobid
	finalize_job
	exit 1
fi

# run any post job scritp(s)
for i in $postjob; do
	[ -f $i ] && runcmd $workdir/$i 
done

# remove the binaries, inputfiles and scripts 
for i in $workrm; do
	[ -f $workdir/$i ] && runcmd rm -rf $workdir/$i
done

# copy what is left back to submit dir
if test $do_copyall = yes; then
	setup_copyall_script
	runcmd $workdir/copyall.sh $workdir $outputdir
else # copy only from the master
	cp -r $workdir/* $outputdir
fi

## remember to cd out of the catalog before you remove it
cd /tmp
runcmd rm -rf $workdir

##  finalize job
finalize_job

#
## The job ends when the script exit.
#
