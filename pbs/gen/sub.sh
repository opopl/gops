#!/bin/bash

source f.sh

#s2 doesn't handle 10, only 6/*{{{*/
#but s2 and s4 can handle 10

# S sub -n 0 -q s2 -nrg 10 chT 2 ; sub -n 1 -q s4 -nrg 10 chT 2

# r -q s1 -od sub_1_chF -radius 2.0 -nrg 1 -nsteps 10000 -temp 0.03 -ap -if 0.0 0.02 0.001 
# 3110-2319
# 381588

#r -q s1 -od sub_0_chF -radius 2.0 -nrg 1 -nsteps 1000 -temp 0.03 -ap -if 0.0 0.02 0.001 
# 3110-2320
# 381589
# f 2324

#r -q s1 -od sub_1_chF -radius 2.0 -nrg 1 -nsteps 10000 -temp 0.03 -ap -if 0.0 0.02 0.001 
# 3110-2326
# 381590

#r -q s1 -od sub_1_chF_nrg-10 -radius 2.0 -nrg 10 -nsteps 10000 -temp 0.03 -ap -if 0.0 0.02 0.001 
# 3110-2327
# 381591

#sub -n 0 -q s2 -nrg 10 chT 3 ; sub -n 1 -q s4 -nrg 10 chT 3
# 3110-2334
# 381592
# 381593

# cpl sub_0_chT_f-0.0 sub_1_chT_f-0.0 sub_0_chF sub_1_chF sub_1_chF_nrg-10
# 0036

# should now define forces + more nsteps
# 1. f_max=0.001

#r -q s1 -od sub_1_chF -radius 2.0 -nrg 1 -nsteps 10000 -temp 0.03 -ap -if 0.0 0.001 0.0001 
# 0111-0049
#381612

#r -q s1 -od sub_2_chF -radius 2.0 -nrg 1 -nsteps 100000 -temp 0.03 -ap -if 0.0 0.001 0.0001 
# 0111-0051
# 381613

# r -q s1  -od sub_3_chF -radius 2.0 -nrg 1 -nsteps 1000000 -temp 0.03 -ap -if 0.0 0.001 0.0001 
# 0111-0111
# 381614

### need to be able to pass walltime as a cmdline parameter..
# 2. f_max= .002

#r -q l1 -od sub_3_chF-2 -radius 2.0 -nrg 1 -nsteps 1000000 -temp 0.03 -ap -if 0.0011 0.002 0.0001 
# 0111-0141
# 381615

#r -q l2 -od sub_3_chF-T-1.0 -radius 2.0 -nrg 1 -nsteps 1000000 -temp 1.0 -ap -if 0.00 0.001 0.0001 
# 0111-0153
# 381616

#r -q l2 -od sub_3_chF-T-0.1 -radius 2.0 -nrg 1 -nsteps 1000000 -temp 0.1 -ap -if 0.00 0.001 0.0001 
# 0111-0217
# 381617

#r -q s2 -od sub_3_chF_nrg-10 -radius 2.0 -nrg 10 -nsteps 1000000 -temp 0.03 -ap -if 0.0 0.001 0.0001 
# 0111-0404
# 381618

#r -q s2 --adj_target 1 -od sub_3_chF_nrg-10_AT-1 -radius 2.0 -nrg 10 -nsteps 1000000 -temp 0.03 -ap -if 0.0 0.003 0.0001 
# 0111-0426
# 381619

# r -q s2 --adj_target 1 -od sub_3_chF_AT-1 -radius 2.0 -nrg 1 -nsteps 1000000 -temp 0.03 -ap -if 0.0 0.003 0.0001 
# 0111-0426
# 381620

# r -q s4 --adj_target 0 -od sub_1_chF_T-1.0-long -radius 2.0 -nrg 1 -nsteps 10000 -temp 1.0 -ap -if 0.0 0.01 0.0001 
# 0111-0445
# 381621

#r -q s4 --adj_target 0 -od sub_1_chF_T-2.0-long -radius 2.0 -nrg 1 -nsteps 10000 -temp 2.0 -ap -if 0.0 0.01 0.0001 
# 0111-0445
# 381622

#r -q s4 --adj_target 0 -od sub_1_chF_T-0.1-long -radius 2.0 -nrg 1 -nsteps 10000 -temp 0.1 -ap -if 0.0 0.01 0.0001 
# 0111-0455
# 381623

# data.old taken in all prev. studies

#SLOPPYCONV 0.0001
#TIGHTCONV  0.0000001
#UPDATES 50 
#DGUESS  10.0
#RADIUS 100.0
#PULL 1 46 0.0
#SAVE 10
#CENTRE
#COMMENT CHANGEACCEPT 50 
#COMMENT CHARMM
#COMMENT DUMP
#COMMENT DEBUG
#COMMENT SORT
#G46
#EDIFF 0.001
#STEPS 9000000 1.0
#MAXBFGS 1.0
#MAXIT 1000 1000
#STEP 1.9 0.0 
#TEMPERATURE 0.035
#TARGET -0.2936531215
#TRACKDATA

#r -q s4 --adj_target 1 -od sub_2_chF_nrg-10_T-0.1_R-5.0_AT-1 \
	#-radius 5.0 -nrg 10 -nsteps 100000 -temp 0.1 -ap -if 0.0 0.03 0.0001 
# 0111-0809
# 381624
# qdel 381624
# 0111-0812
# 381625

#r -q s4 --adj_target 0 -od sub_2_chF_nrg-10_T-0.1_R-5.0 \
	#-radius 5.0 -nrg 10 -nsteps 100000 -temp 0.1 -ap -if 0.0 0.03 0.0001 
# 0111-0812
# 381626

#sub -n 0 -q s2 -nrg 10 chT 4
#sub -n 1 -q s4 -nrg 10 chT 4
# 0111-0820
#381627
#381628
#381629
#381630
#381631
#381632
#381633
#381634
#381635
#381636
#381637
#381638
#381639
#381640
#381641
#381642
#381643
#381644
#381645
#381646

# 0111 - cancelled all my remaining jobs.
# concentrate on n=10000, and very small refinement for forces.

# r -q s1 --adj_target 0 -od sub_1_chF -radius 5.0 -nrg 10 -nsteps 10000 -temp 0.01 -ap -if 0.0 0.01 0.0001
# 0111-1116
# 381647 

#r -q s1 --adj_target 0 -od sub_0_chF -radius 5.0 -nrg 10 -nsteps 1000 -temp 0.01 -ap -if 0.0 0.01 0.0001
# 0111-1123
# 381648 

#r -q s2 --adj_target 0 -od sub_2_chF_nrg-1 -radius 5.0 -nsteps 100000 -temp 0.01 -ap -if 0.0 0.01 0.0001
# 0111-1138
# 381649 

#r -q s2 --adj_target 0 -od sub_2_chF_nrg-10 -radius 5.0 -nsteps 100000 -temp 0.01 -ap -if 0.0 0.01 0.0001
# 0111-1138
# 381650 

#
#r -q s2 --adj_target 0 -od sub_1_chF_nrg-1 -radius 5.0 -nsteps 10000 -temp 0.01 -ap -if 0.0 0.01 0.0001
#r -q s2 --adj_target 0 -od sub_1_chF_nrg-10 -radius 5.0 -nsteps 10000 -temp 0.01 -ap -if 0.0 0.01 0.0001
#r -q s4 --adj_target 0 -od sub_2_chF_nrg-1 -radius 5.0 -nsteps 100000 -temp 0.01 -ap -if 0.0 0.01 0.0001
#r -q s4 --adj_target 0 -od sub_2_chF_nrg-10 -radius 5.0 -nsteps 100000 -temp 0.01 -ap -if 0.0 0.01 0.0001
# 0111-1159
#381651.clust.ch.cam.ac.uk
#381652.clust.ch.cam.ac.uk
#381653.clust.ch.cam.ac.uk
#381654.clust.ch.cam.ac.uk
#
# larger T=0.05

#r -q s2 --adj_target 0 -od sub_1_chF_nrg-1_T-0.05 -radius 5.0 -nsteps 10000 -temp 0.05 -ap -if 0.0 0.01 0.0001
#r -q s2 --adj_target 0 -od sub_1_chF_nrg-10_T-0.05 -radius 5.0 -nsteps 10000 -temp 0.05 -ap -if 0.0 0.01 0.0001
#r -q s4 --adj_target 0 -od sub_2_chF_nrg-1_T-0.05 -radius 5.0 -nsteps 100000 -temp 0.05 -ap -if 0.0 0.01 0.0001
#r -q s4 --adj_target 0 -od sub_2_chF_nrg-10_T-0.05 -radius 5.0 -nsteps 100000 -temp 0.05 -ap -if 0.0 0.01 0.0001

#381655.clust.ch.cam.ac.uk
#381656.clust.ch.cam.ac.uk
#381657.clust.ch.cam.ac.uk
#381658.clust.ch.cam.ac.uk

# at office now, something strange happened with the directories,
# have to re-submit now
# 0111-1555

#
# don't forget about the correspondence:
# 
# 	s1 =>	1000
# 	s2 =>	10000
# 	s4 =>	100000
#

n_sleep=5

#temp=0.05
#r -q s2 --adj_target 0 -od sub_1_chF_nrg-1_T-0.05 -radius 5.0 -nsteps 10000 -temp 0.05 -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#r -q s2 --adj_target 0 -od sub_1_chF_nrg-10_T-0.05 -radius 5.0 -nsteps 10000 -temp 0.05 -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#r -q s4 --adj_target 0 -od sub_2_chF_nrg-1_T-0.05 -radius 5.0 -nsteps 100000 -temp 0.05 -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#r -q s4 --adj_target 0 -od sub_2_chF_nrg-10_T-0.05 -radius 5.0 -nsteps 100000 -temp 0.05 -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#
#381659.clust.ch.cam.ac.uk
#381660.clust.ch.cam.ac.uk
#381661.clust.ch.cam.ac.uk
#381662.clust.ch.cam.ac.uk
# deleted, added n_sleep, resubmit now

# 1604
#381663.clust.ch.cam.ac.uk
#381664.clust.ch.cam.ac.uk
#381665.clust.ch.cam.ac.uk
#381666.clust.ch.cam.ac.uk

#temp=0.01
#r -q s2 --adj_target 0 -od sub_1_chF_nrg-1_T-$temp -radius 5.0 -nsteps 10000 -temp $temp -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#r -q s2 --adj_target 0 -od sub_1_chF_nrg-10_T-$temp -radius 5.0 -nsteps 10000 -temp $temp -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#r -q s4 --adj_target 0 -od sub_2_chF_nrg-1_T-$temp -radius 5.0 -nsteps 100000 -temp $temp -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#r -q s4 --adj_target 0 -od sub_2_chF_nrg-10_T-$temp -radius 5.0 -nsteps 100000 -temp $temp -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#381668.clust.ch.cam.ac.uk
#381669.clust.ch.cam.ac.uk
#381670.clust.ch.cam.ac.uk
#381671.clust.ch.cam.ac.uk

#temp=1.0
#r -q s2 --adj_target 0 -od sub_1_chF_nrg-1_T-$temp -radius 5.0 -nsteps 10000 -temp $temp -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#r -q s2 --adj_target 0 -od sub_1_chF_nrg-10_T-$temp -radius 5.0 -nsteps 10000 -temp $temp -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#r -q s4 --adj_target 0 -od sub_2_chF_nrg-1_T-$temp -radius 5.0 -nsteps 100000 -temp $temp -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#r -q s4 --adj_target 0 -od sub_2_chF_nrg-10_T-$temp -radius 5.0 -nsteps 100000 -temp $temp -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#381673.clust.ch.cam.ac.uk
#381674.clust.ch.cam.ac.uk
#381675.clust.ch.cam.ac.uk
#381676.clust.ch.cam.ac.uk
# submitted, time: 1615

# with ORGYR keyword

#temp=0.01
#r -q s4 --adj_target 0 -od sub_2_chF_nrg-1_T-"$temp"_R-2.0-ORGYR -radius 2.0 -nsteps 100000 -temp $temp -ap -if 0.0 0.01 0.0001
#sleep $n_sleep

#
# submitted 1640
# 381677
#

#temp=0.01
#r -q s4 --adj_target 0 -od sub_2_chF_nrg-1_T-"$temp"_R-2.0-ORGYR -radius 2.0 -nsteps 100000 -temp $temp -ap -if 0.0 0.01 0.0001
#sleep $n_sleep/*}}}*/

n_sleep=5
temp=0.05
nrg=1

# at office now: 1853

#r -q s2 --adj_target 0 -nrg $nrg -od sub_1_chF_NRG-1_T-$temp -radius 5.0 -nsteps 10000 -temp $temp -ap -if 0.0 0.01 0.0001

#sleep $n_sleep

# 381679
#temp=0.01
#r -q s4 --adj_target 0 -nrg $nrg -od sub_1_chF_NRG-1_T-$temp -radius 5.0 -nsteps 10000 -temp $temp -ap -if 0.0 0.01 0.0001

# 381680

#sub_1_chF_NRG-1_T-0.01:
# 0111-1908
# 0.0003 - seems wrong. Refine...
#temp=0.01
#nrg=1
#force=0.0003
#queue=s4
#r -q $queue -nrg $nrg -od f-3E-4 -radius 5.0 -nsteps 10000 -temp 0.01 -ap -force $force 
#r -q $queue -nrg $nrg -od f-3E-4 -radius 5.0 -nsteps 10000 -temp 0.05 -ap -force $force 
#r -q $queue -nrg $nrg -od f-3E-4 -radius 5.0 -nsteps 50000 -temp 0.01 -ap -force $force 

# 1914, 3 jobs: 3E-4
#381681.clust.ch.cam.ac.uk
#381682.clust.ch.cam.ac.uk
#381683.clust.ch.cam.ac.uk
# resubmit, was -d instead of -od...
# reset queue to s1...
#queue=s1
#od=f-3E-4
#r -q $queue -nrg $nrg -od $od -radius 5.0 -nsteps 10000 -temp 0.01 -ap -force $force 
#sleep $n_sleep
#r -q $queue -nrg $nrg -od $od -radius 5.0 -nsteps 10000 -temp 0.05 -ap -force $force 
#sleep $n_sleep
#r -q $queue -nrg $nrg -od $od -radius 5.0 -nsteps 50000 -temp 0.01 -ap -force $force 
#sleep $n_sleep
#r -q $queue -nrg $nrg -od $od -radius 5.0 -nsteps 1000000 -temp 0.01 -ap -force $force 
#sleep $n_sleep
# seems there is a phase structural transition(s) ! 
# ...  between 2E-4 and 3E-4 ?
# 381688, ..., 381691  
# resubmit, was dn instead of dnj in loop_forces.sh
# 0111-1927: resubmitting...
# 381691, ... , 95 
#
# resubmit...
# ... 96-99
# seems finally is running
# resubmit again, to clean all the output files...
# 19-34, resubmitting...
# 381700, ..., 381703
# done

#
# now, we should investigate the region where we suspect the phase structural transition
# 0111-1948
# between 2E-4 and 3E-4, in 1E-5 increments, n=10000 t=0.01 should be enough?

#od=f-3E-4-left
#queue=s8
#nrg=1
#r -q $queue -nrg $nrg -od $od -radius 5.0 -nsteps 10000 -temp 0.01 -ap -if 0.0003 0.0004 0.00001 
#sleep $n_sleep
#
# 381704 
# 1952 submitted
# again...

# od=f-3E-4-left
# queue=s8
# nrg=1
# r -q $queue -nrg $nrg -od $od -radius 5.0 -nsteps 10000 -temp 0.01 -ap -if 0.0002 0.0003 0.00001 
# sleep $n_sleep
#
# 381705 
# 1955 submitted
#

od="sub_2_chF" 
nrg=1
temp=0.03
nsteps=100000
queue=s8
radius=5.0
n_sleep=3

#force_min=0.0 ; force_max=0.001 ; force_i=0.0001 ; source define_forces_0.sh

forces=( 0.001 0.005 0.01 0.05 0.1 ) 

let i=0
for force in "${forces[@]}"
	do
		r -q $queue -nrg $nrg -od $od -radius $radius -nsteps $nsteps -temp $temp -ap -force $force 
		sleep $n_sleep
		i=$(($i+1))
done

# submitted to s1, 0111-2249
# 381706, ... , 710
# doesn't want to run !!!
#
