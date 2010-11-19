
# adjust the GMIN input file 'data'

change.sh pull_start $pull_start
change.sh pull_end $pull_end
change.sh nsteps $nsteps
change.sh force $force
change.sh temperature $temp
change.sh radius $radius
change.sh sys $sys
echo "COMMENT Changed by O.P. on `date_dm_hm` " >> data.G46
